#include <iostream>
#include <fstream>
#include <vector>
#include <functional>

#include "bitstream.h"
#include "tables.h"
#include "wav.h"

#define NOT_IMPLEMENTED() do { throw std::exception("not implemented"); } while(false);

enum class Mode { Stereo, IntensityStereo, DualChannel, SingleChannel };

enum class Emphasis { None, Microseconds_50_15, Reserved, CCITT_J17 };

struct Header {
    int id;
    int layer;
    bool protection_bit;
    int bitrate;
    int sampling_frequency;
    bool padding_bit;
    bool private_bit;
    Mode mode;
    int mode_extension;
    bool copyright_bit;
    bool original_bit;
    Emphasis emphasis;
    unsigned short crc16;
    int frame_size;
};

struct AudioDataI {
    int allocations[2][32] = {};
    int levels[2][32] = {};
    int scale_factors[2][32] = {};
    int samples[2][32][12] = {};
    float requantized_samples[2][32][12] = {};
};

struct AudioDataII {
    int allocations[2][32] = {};
    int scfsi[2][32] = {};
    float scale_factors[2][32][3] = {};
    int samples[2][32][36] = {};
    float requantized_samples[2][32][36] = {};
};

struct AudioDataIII {
    int samples[2][576] = {};
};

float V[1024] = {};

size_t left_most_bit_index(int value) {
    size_t counter = 0;
    while (value > 0) {
        value >>= 1;
        ++counter;
    }
    return counter;
}

bool synchronize(BitStream& buffer) {
    size_t one_counter = 0;
    while (one_counter < 12) {
        if (buffer.eof()) return false;
        uint8_t bit = buffer.read_bit();
        one_counter = bit ? one_counter + 1 : 0;
        if (!bit) buffer.go_to_next_byte();
    }
    return true;
}

float requantize_II(int raw_sample, int num_bits, int quant_index) {
    const int max_value = 1 << (num_bits - 1);
    const int bit_mask = max_value - 1;

    float fraction = (raw_sample & bit_mask) / static_cast<float>(max_value);
    if (!(raw_sample & max_value)) fraction += -1;

    return layer_II_quantization_class_C[quant_index] * (fraction + layer_II_quantization_class_D[quant_index]);
}

float requantize_I(int raw_sample, int num_bits) {
    const int max_value = 1 << (num_bits - 1);
    const int bit_mask = max_value - 1;
    const float residue = 1.0f / max_value;

    float fraction = (raw_sample & bit_mask) / static_cast<float>(max_value) + residue;
    if (!(raw_sample & max_value)) fraction += -1;

    const float a = 1 << num_bits;
    const float b = (1 << num_bits) - 1;
    const float scalar = a / b;

    return scalar * fraction;
}

Header read_header(BitStream& bitstream) {
    Header header;
    header.id                   = bitstream.read_bit();
    header.layer                = table_layer[bitstream.read_bits(2)];
    header.protection_bit       = bitstream.read_bit();
    header.bitrate              = table_bitrate_per_layer[header.layer - 1][bitstream.read_bits(4)];
    header.sampling_frequency   = table_samplerate[bitstream.read_bits(2)];
    header.padding_bit          = bitstream.read_bit();
    header.private_bit          = bitstream.read_bit();
    header.mode                 = static_cast<Mode>(bitstream.read_bits(2));
    header.mode_extension       = bitstream.read_bits(2);
    header.copyright_bit        = bitstream.read_bit();
    header.original_bit         = bitstream.read_bit();
    header.emphasis             = static_cast<Emphasis>(bitstream.read_bits(2));
    if (!header.protection_bit)
        header.crc16            = bitstream.read_bits<unsigned short>(16);
    header.frame_size           = 144 * header.bitrate * 1000 / header.sampling_frequency + header.padding_bit;
    return header;
}

struct layer_II_quantization_table_info {
    int (*table)[32][16];
    int (*nbal)[32];
    int sblimit;
};

layer_II_quantization_table_info layer_II_get_quantization_table(const Header& header) {
    switch (header.sampling_frequency) {
    case 48000:
        switch (header.bitrate) {
        case 56: case 64: case 80: case 96: case 112: case 128: case 160: case 192:
            return { &layer_II_quantizations_a, &layer_II_quantizations_a_nbal, layer_II_quantizations_a_sblimit };
        case 32: case 48:
            return { &layer_II_quantizations_c, &layer_II_quantizations_c_nbal, layer_II_quantizations_c_sblimit };
        }
        break;
    case 44100:
        switch (header.bitrate) {
        case 56: case 64: case 80:
            return { &layer_II_quantizations_a, &layer_II_quantizations_a_nbal, layer_II_quantizations_a_sblimit };
        case 96: case 112: case 128: case 160: case 192:
            return { &layer_II_quantizations_b, &layer_II_quantizations_b_nbal, layer_II_quantizations_b_sblimit };
        case 32: case 48:
            return { &layer_II_quantizations_c, &layer_II_quantizations_c_nbal, layer_II_quantizations_c_sblimit };
        }
    case 32000:
        switch (header.bitrate) {
        case 56: case 64: case 80:
            return { &layer_II_quantizations_a, &layer_II_quantizations_a_nbal, layer_II_quantizations_a_sblimit };
        case 96: case 112: case 128: case 160: case 192:
            return { &layer_II_quantizations_b, &layer_II_quantizations_b_nbal, layer_II_quantizations_b_sblimit };
        case 32: case 48:
            return { &layer_II_quantizations_d, &layer_II_quantizations_d_nbal, layer_II_quantizations_d_sblimit };
        }
    default:
        break;
    }

    throw std::exception();
}

void loop_over_sb_and_ch(int channels, int bound, int bound_limit, std::function<void(int ch, int sb, bool is_intensity)> f) {
    for (size_t sb = 0; sb < bound; sb++) {
        for (size_t ch = 0; ch < channels; ch++) {
            f(ch, sb, false);
        }
    }

    for (size_t sb = bound; sb < bound_limit; sb++) {
        f(0, sb, true);
    }
}

const HuffmanEntry4& huffman_decode(RingBitStream& bitstream, const HuffmanEntry4* table, int size) {
    int result = bitstream.read_bit();
    int bits_read = 1;

    while (true) {
        for (size_t i = 0; i < size; i++) {
            if (table[i].hcod == result && table[i].hlen == bits_read) {
                return table[i];
            }
        }

        result = (result << 1) | bitstream.read_bit();
        ++bits_read;
    }

    throw std::exception();
}

const HuffmanEntry2& huffman_decode(RingBitStream& bitstream, const HuffmanEntry2* table, int size) {
    if (size == 1) {
        return *table;
    }

    int bits[22] = {};
    int result = bitstream.read_bit();
    int bits_read = 1;
    bits[0] = result;

    while (bits_read < 22) {
        for (size_t i = 0; i < size; i++) {
            if (table[i].hcod == result && table[i].hlen == bits_read) {
                return table[i];
            }
        }

        result = (result << 1) | bitstream.read_bit();
        bits[bits_read] = result & 1;
        ++bits_read;
    }

    throw std::exception();
}

AudioDataIII read_audio_data_III(const Header& header, BitStream& bitstream, RingBitStream& reservoir) {
    int private_bits = {};
    int scfsi[2][4] = {};
    int part2_3_length[2][2] = {};
    int big_values[2][2] = {};
    int global_gain[2][2] = {};
    int scalefac_compress[2][2] = {};
    int window_switching_flag[2][2] = {};
    int block_type[2][2] = {};
    int mixed_block_flag[2][2] = {};
    int table_select[2][2][3] = {};
    int sub_block_gain[2][2][3] = {};
    int region0_count[2][2] = {};
    int region1_count[2][2] = {};
    int preflag[2][2] = {};
    int scalefac_scale[2][2] = {};
    int count1table_select[2][2] = {};

    reservoir.seek_to_end();

    const int channels = header.mode == Mode::SingleChannel ? 1 : 2;
    const int main_data_begin = bitstream.read_bits(9);

    if (header.mode == Mode::SingleChannel) {
        private_bits = bitstream.read_bits(5);
    } else {
        private_bits = bitstream.read_bits(3);
    }

    for (size_t ch = 0; ch < channels; ch++) {
        for (size_t scfsi_band = 0; scfsi_band < 4; scfsi_band++) {
            scfsi[ch][scfsi_band] = bitstream.read_bits(1);
        }
    }

    for (size_t gr = 0; gr < 2; gr++) {
        for (size_t ch = 0; ch < channels; ch++) {
            part2_3_length[gr][ch] = bitstream.read_bits(12);
            big_values[gr][ch] = bitstream.read_bits(9);
            global_gain[gr][ch] = bitstream.read_bits(8);
            scalefac_compress[gr][ch] = bitstream.read_bits(4);
            window_switching_flag[gr][ch] = bitstream.read_bits(1);
            if (window_switching_flag[gr][ch]) {
                block_type[gr][ch] = bitstream.read_bits(2);
                mixed_block_flag[gr][ch] = bitstream.read_bits(1);
                for (size_t region = 0; region < 2; region++)
                    table_select[gr][ch][region] = bitstream.read_bits(5);
                for (size_t window = 0; window < 3; window++)
                    sub_block_gain[gr][ch][window] = bitstream.read_bits(3);
                region0_count[gr][ch] = (block_type[gr][ch] == 2 && !mixed_block_flag[gr][ch]) ? 8 : 7;
                region1_count[gr][ch] = (block_type[gr][ch] == 2 && !mixed_block_flag[gr][ch]) ? 36 : (20 - region0_count[gr][ch]);
            } else {
                for (size_t region = 0; region < 3; region++)
                    table_select[gr][ch][region] = bitstream.read_bits(5);
                region0_count[gr][ch] = bitstream.read_bits(4);
                region1_count[gr][ch] = bitstream.read_bits(3);
            }
            preflag[gr][ch] = bitstream.read_bits(1);
            scalefac_scale[gr][ch] = bitstream.read_bits(1);
            count1table_select[gr][ch] = bitstream.read_bits(1);
        }
    }

    _ASSERT(bitstream.get_current_bit() == 0);

    for (size_t i = 17 + 4; i < header.frame_size; i++) {
        unsigned char data = bitstream.read_bits<unsigned char>(8);
        reservoir.write(&data, 1);
    }

    if (main_data_begin > 0) {
        reservoir.seek_relative(-main_data_begin);
    }

    int scalefac_bits_count[2] = { 0 };

    int scalefac_l[2][2][21] = {};
    int scalefac_s[2][2][8][3] = {};

    for (size_t gr = 0; gr < 2; gr++) {
        int current_position = reservoir.position();
        for (size_t ch = 0; ch < channels; ch++) {
            if (window_switching_flag[gr][ch] == 1 && block_type[gr][ch] == 2) {
                if (mixed_block_flag[gr][ch]) {
                    for (size_t sfb = 0; sfb < 8; sfb++) {
                        const int bits = layer_III_scalefac_compress_slen1[scalefac_compress[gr][ch]];
                        scalefac_l[gr][ch][sfb] = reservoir.read_bits(bits);
                    }
                    for (size_t sfb = 3; sfb < 12; sfb++) {
                        for (size_t window = 0; window < 3; window++) {
                            const int bits = sfb <= 5 ? layer_III_scalefac_compress_slen1[scalefac_compress[gr][ch]] :
                                                        layer_III_scalefac_compress_slen2[scalefac_compress[gr][ch]];
                            scalefac_s[gr][ch][sfb][window] = reservoir.read_bits(bits);
                        }
                    }
                } else {
                    for (size_t sfb = 0; sfb < 12; sfb++) {
                        for (size_t window = 0; window < 3; window++) {
                            const int bits = sfb <= 5 ? layer_III_scalefac_compress_slen1[scalefac_compress[gr][ch]] :
                                                        layer_III_scalefac_compress_slen2[scalefac_compress[gr][ch]];
                            scalefac_s[gr][ch][sfb][window] = reservoir.read_bits(bits);
                        }
                    }
                }
            } else {
                auto scalefactor_bits = [&gr, &ch, &scalefac_compress](int sfb) -> int {
                    return sfb <= 10 ? layer_III_scalefac_compress_slen1[scalefac_compress[gr][ch]]
                                     : layer_III_scalefac_compress_slen2[scalefac_compress[gr][ch]];
                };

                if ((scfsi[ch][0] == 0) || (gr == 0)) {
                    for (size_t sfb = 0; sfb < 6; sfb++) {
                        const int bits = scalefactor_bits(sfb);
                        scalefac_l[gr][ch][sfb] = reservoir.read_bits(bits);
                    }
                }
                if ((scfsi[ch][1] == 0) || (gr == 0)) {
                    for (size_t sfb = 6; sfb < 11; sfb++) {
                        const int bits = scalefactor_bits(sfb);
                        scalefac_l[gr][ch][sfb] = reservoir.read_bits(bits);
                    }
                }
                if ((scfsi[ch][2] == 0) || (gr == 0)) {
                    for (size_t sfb = 11; sfb < 16; sfb++) {
                        const int bits = scalefactor_bits(sfb);
                        scalefac_l[gr][ch][sfb] = reservoir.read_bits(bits);
                    }
                }
                if ((scfsi[ch][3] == 0) || (gr == 0)) {
                    for (size_t sfb = 16; sfb < 21; sfb++) {
                        const int bits = scalefactor_bits(sfb);
                        scalefac_l[gr][ch][sfb] = reservoir.read_bits(bits);
                    }
                }
            }
        }
        scalefac_bits_count[gr] = reservoir.position() - current_position;
    }

    AudioDataIII result;

    for (size_t gr = 0; gr < 2; gr++) {
        int count = 0;

        for (size_t ch = 0; ch < channels; ch++) {
            int scaleFactorBandIndex1 = region0_count[gr][ch] + 1;
            int scaleFactorBandIndex2 = scaleFactorBandIndex1 + region1_count[gr][ch] + 1;
            if (scaleFactorBandIndex2 >= ScaleFactorBandsLong[header.sampling_frequency].size()) {
                scaleFactorBandIndex2 = ScaleFactorBandsLong[header.sampling_frequency].size() - 1;
            }

            int region1_start = ScaleFactorBandsLong[header.sampling_frequency][scaleFactorBandIndex1].start;
            int region2_start = ScaleFactorBandsLong[header.sampling_frequency][scaleFactorBandIndex2].start;
            if (window_switching_flag[gr][ch] && block_type[gr][ch] == 2) {
                region1_start = 36;
                region2_start = 576;
            }

            int current_position = reservoir.position();

            const HuffmanTable2* table = nullptr;
            size_t i = 0;
            for (; i < big_values[gr][ch] * 2; i += 2) {
                if (i < region1_start) {
                    table = &HuffmanTables2[table_select[gr][ch][0]];
                } else if (i < region2_start) {
                    table = &HuffmanTables2[table_select[gr][ch][1]];
                } else {
                    table = &HuffmanTables2[table_select[gr][ch][2]];
                }

                if (table == nullptr || table->table == nullptr) {
                    std::cout << "UNKNOWN TABLE: " << table_select[gr][ch][0] << std::endl;
                    throw std::exception();
                }

                const HuffmanEntry2& huff = huffman_decode(reservoir, table->table, table->table_size);
                int linbitsx = huff.x;
                int linbitsy = huff.y;
                if (linbitsx == 15 && table->linbits > 0) {
                    linbitsx += reservoir.read_bits(table->linbits);
                }
                if (linbitsx != 0) {
                    bool sign = reservoir.read_bit();
                    if (sign) linbitsx = -linbitsx;
                }
                if (linbitsy == 15 && table->linbits > 0) {
                    linbitsy += reservoir.read_bits(table->linbits);
                }
                if (linbitsy != 0) {
                    bool sign = reservoir.read_bit();
                    if (sign) linbitsy = -linbitsy;
                }

                result.samples[gr][count++] = linbitsx;
                result.samples[gr][count++] = linbitsy;
            }

            const HuffmanEntry4* count1table = count1table_select[gr][ch] ? huffman_table_B : huffman_table_A;
            int granule_bits_read = (reservoir.position() - current_position) + scalefac_bits_count[gr];

            // count1 is not known. We have to read huffman encoded values
            // until we've exhausted the granule's bits. We know the size of
            // the granule from part2_3_length, which is the number of bits
            // used for scaleactors and huffman data (in the granule).
            while (granule_bits_read < part2_3_length[gr][ch] && count < 576) {
                const HuffmanEntry4& entry = huffman_decode(reservoir, count1table, 16);
                int v = entry.v;
                if (v != 0 && reservoir.read_bit()) v = -v;
                int w = entry.w;
                if (w != 0 && reservoir.read_bit()) w = -w;
                int x = entry.x;
                if (x != 0 && reservoir.read_bit()) x = -x;
                int y = entry.y;
                if (y != 0 && reservoir.read_bit()) y = -y;

                result.samples[gr][count++] = v;
                result.samples[gr][count++] = w;
                result.samples[gr][count++] = x;
                result.samples[gr][count++] = y;

                granule_bits_read = (reservoir.position() - current_position) + scalefac_bits_count[gr];
            }

            if (granule_bits_read > part2_3_length[gr][ch]) {
                reservoir.rewind(granule_bits_read - part2_3_length[gr][ch]);
                count--;
            }

            if (granule_bits_read < part2_3_length[gr][ch]) {
                for (size_t i = granule_bits_read; i < part2_3_length[gr][ch]; i++) {
                    reservoir.read_bit();
                }
            }
        }
    }

    return result;
}

AudioDataII read_audio_data_II(const Header& header, BitStream& bitstream) {
    AudioDataII data;

    auto quantization_info = layer_II_get_quantization_table(header);

    const int channels = header.mode == Mode::SingleChannel ? 1 : 2;
    const int sblimit = quantization_info.sblimit;
    const int bound = sblimit; // TODO: for Joint Stereo this is not true

    loop_over_sb_and_ch(channels, bound, sblimit, [&](int ch, int sb, bool is_intensity) {
        const int al = bitstream.read_bits((*quantization_info.nbal)[sb]);
        data.allocations[ch][sb] = (*quantization_info.table)[sb][al];
        if (is_intensity) data.allocations[1][sb] = data.allocations[0][sb];
    });

    loop_over_sb_and_ch(channels, sblimit, sblimit, [&](int ch, int sb, bool is_intensity) {
        if (data.allocations[ch][sb] == 0) return;
        data.scfsi[ch][sb] = bitstream.read_bits(2);
    });

    loop_over_sb_and_ch(channels, sblimit, sblimit, [&](int ch, int sb, bool is_intensity) {
        if (data.allocations[ch][sb] == 0) return;
        switch (data.scfsi[ch][sb]) {
        case 0:
            data.scale_factors[ch][sb][0] = table_scale_factors[bitstream.read_bits(6)];
            data.scale_factors[ch][sb][1] = table_scale_factors[bitstream.read_bits(6)];
            data.scale_factors[ch][sb][2] = table_scale_factors[bitstream.read_bits(6)];
            break;
        case 2:
            data.scale_factors[ch][sb][0] = table_scale_factors[bitstream.read_bits(6)];
            data.scale_factors[ch][sb][1] = data.scale_factors[ch][sb][0];
            data.scale_factors[ch][sb][2] = data.scale_factors[ch][sb][0];
            break;
        case 1:
            data.scale_factors[ch][sb][0] = table_scale_factors[bitstream.read_bits(6)];
            data.scale_factors[ch][sb][1] = data.scale_factors[ch][sb][0];
            data.scale_factors[ch][sb][2] = table_scale_factors[bitstream.read_bits(6)];
            break;
        case 3:
            data.scale_factors[ch][sb][0] = table_scale_factors[bitstream.read_bits(6)];
            data.scale_factors[ch][sb][1] = table_scale_factors[bitstream.read_bits(6)];
            data.scale_factors[ch][sb][2] = data.scale_factors[ch][sb][1];
            break;
        }
    });

    // TODO: Theres probably a smarter way to do this
    auto quantization_class_index = [](int value) -> int {
        for (size_t i = 0; i < 17; i++) {
            if (layer_II_quantization_class_num_steps[i] == value) {
                return i;
            }
        }
        throw std::exception();
    };

    for (size_t gr = 0; gr < 12; gr++) {
        loop_over_sb_and_ch(channels, bound, sblimit, [&](int ch, int sb, bool is_intensity) {
            if (data.allocations[ch][sb] == 0) return;

            const int quant_index = quantization_class_index(data.allocations[ch][sb]);
            const int code_width = layer_II_quantization_class_bits_per_cw[quant_index];

            if (layer_II_quantization_class_group[quant_index]) {
                const int nlevels = data.allocations[ch][sb];
                const size_t sample_code_highest_bit_index = left_most_bit_index(nlevels);

                int sample_code = bitstream.read_bits(code_width);

                for (size_t i = 0; i < 3; i++) {
                    data.samples[ch][sb][3 * gr + i] = sample_code % nlevels;
                    sample_code = sample_code / nlevels;

                    const float scale_factor = data.scale_factors[ch][sb][gr / 4];
                    const int sample = data.samples[ch][sb][3 * gr + i];
                    const float requantized = requantize_II(sample, sample_code_highest_bit_index, quant_index);

                    data.requantized_samples[ch][sb][3 * gr + i] = scale_factor * requantized;
                    if (is_intensity) data.requantized_samples[1][sb][3 * gr + i] = data.requantized_samples[0][sb][3 * gr + i];
                }
            } else {
                for (size_t s = 0; s < 3; s++) {
                    data.samples[ch][sb][3 * gr + s] = bitstream.read_bits(code_width);

                    const float scale_factor = data.scale_factors[ch][sb][gr / 4];
                    const int sample = data.samples[ch][sb][3 * gr + s];
                    const float requantized = requantize_II(sample, code_width, quant_index);

                    data.requantized_samples[ch][sb][3 * gr + s] = scale_factor * requantized;
                    if (is_intensity) data.requantized_samples[1][sb][3 * gr + s] = data.requantized_samples[0][sb][3 * gr + s];
                }
            }
        });
    }

    return data;
}

AudioDataI read_audio_data_I(const Header& header, BitStream& bitstream) {
    AudioDataI data;

    const int channels = header.mode == Mode::SingleChannel ? 1 : 2;
    const int bound = header.mode == Mode::IntensityStereo ? table_bounds[header.mode_extension] : 32;

    loop_over_sb_and_ch(channels, bound, 32, [&](int ch, int sb, bool is_intensity) {
        data.allocations[ch][sb] = table_allocation[bitstream.read_bits(4)];
        if (is_intensity) data.allocations[1][sb] = data.allocations[0][sb];
    });

    loop_over_sb_and_ch(channels, 32, 32, [&](int ch, int sb, bool is_intensity) {
        if (data.allocations[ch][sb] == 0) return;
        data.scale_factors[ch][sb] = bitstream.read_bits(6);
    });

    for (size_t s = 0; s < 12; s++) {
        loop_over_sb_and_ch(channels, bound, 32, [&](int ch, int sb, bool is_intensity) {
            if (data.allocations[ch][sb] == 0) return;
            const int num_bits = data.allocations[ch][sb];
            const float scale_factor = table_scale_factors[data.scale_factors[ch][sb]];
            data.samples[ch][sb][s] = bitstream.read_bits(num_bits);
            data.requantized_samples[ch][sb][s] = scale_factor * requantize_I(data.samples[ch][sb][s], num_bits);
            if (is_intensity) data.requantized_samples[1][sb][s] = data.requantized_samples[0][sb][s];
        });
    }

    return data;
}

// ISO/IEC 11172-3 (Figure A.2)
void synthesis(float samples[32], float result[32]) {
    for (size_t i = 1023; i >= 64; i--) {
        V[i] = V[i - 64];
    }

    for (size_t i = 0; i < 64; i++) {
        V[i] = 0;
        for (size_t k = 0; k < 32; k++) {
            const float N = cos((16 + i) * (2 * k + 1) * M_PI / 64.0f);
            V[i] += N * samples[k];
        }
    }

    float U[512];
    for (size_t i = 0; i < 8; i++) {
        for (size_t j = 0; j < 32; j++) {
            U[i * 64 + j] = V[i * 128 + j];
            U[i * 64 + 32 + j] = V[i * 128 + 96 + j];
        }
    }

    float W[512];
    for (size_t i = 0; i < 512; i++) {
        W[i] = U[i] * table_window[i];
    }

    for (size_t j = 0; j < 32; j++) {
        result[j] = 0;
        for (size_t k = 0; k < 16; k++) {
            result[j] += W[j + 32 * k];
        }
    }
}

short clamp(float value) {
    if (value > 32767) return 32767;
    if (value < -32768) return -32768;
    return value;
}

int main()
{
    char* content;
    size_t content_length;
    {
        std::ifstream infile(R"(..\res\chirp.mp3)", std::ios::binary);
        infile.seekg(0, std::ios::end);
        content_length = infile.tellg();
        infile.seekg(0, std::ios::beg);

        content = new char[content_length];
        infile.read(content, content_length);
    }

    auto bitstream = BitStream(content, content_length);

    std::ofstream outfile;
    outfile.open(R"(..\res\out.wav)", std::ios::binary);
    outfile.write(new char[44], 44);
    size_t sample_count = 0;

    float in_samples[32];
    float out_samples[32];

    auto write_samples = [&out_samples, &outfile, &sample_count]() {
        for (size_t j = 0; j < 32; j++) {
            const short sample = clamp(32767 * out_samples[j]);
            outfile.write(reinterpret_cast<const char*>(&sample), 2);
        }
        sample_count += 32;
    };


    RingBitStream reservoir { 65536 };

    std::ofstream f;
    f.open("iii.csv");

    while (!bitstream.eof()) {
        if (!synchronize(bitstream)) {
            std::cout << "no sync found!" << std::endl;
            break;
        }
        std::cout << "Found frame @ " << bitstream.get_current_byte() << std::endl;

        const size_t frame_position = bitstream.get_current_byte() - 1;

        Header header = read_header(bitstream);
        std::cout << "Layer " << header.layer << ", " << header.bitrate << ", " << header.padding_bit << ", " << header.protection_bit << ", " << header.frame_size << std::endl;
        if (header.layer == 1) {
            AudioDataI audioData = read_audio_data_I(header, bitstream);

            for (size_t i = 0; i < 12; i++) {
                for (size_t j = 0; j < 32; j++) {
                    in_samples[j] = audioData.requantized_samples[0][j][i];
                }
                synthesis(in_samples, out_samples);
                write_samples();
            }
        } else if (header.layer == 2) {
            AudioDataII audioData = read_audio_data_II(header, bitstream);

            for (size_t i = 0; i < 36; i++) {
                for (size_t j = 0; j < 32; j++) {
                    in_samples[j] = audioData.requantized_samples[0][j][i];
                }
                synthesis(in_samples, out_samples);
                write_samples();
            }
        } else if (header.layer == 3) {
            AudioDataIII audioData = read_audio_data_III(header, bitstream, reservoir);
            for (int gr = 0; gr < 2; gr++) {
                for (int i = 0; i < 576; i++) {
                    f << audioData.samples[gr][i];
                    if (i < 575) f << ",";
                }
                f << "\n";
            }
            //break;
        }
        else {
            bitstream.go_to_byte(bitstream.get_current_byte() - 2);
        }
    }

    write_wav_header(outfile, 44100, sample_count);
}
