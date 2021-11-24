#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <functional>

#include "bitstream.h"
#include "tables.h"
#include "wav.h"

#define NOT_IMPLEMENTED() do { throw std::exception("not implemented"); } while(false);

enum class Mode { Stereo, IntensityStereo, DualChannel, SingleChannel };

enum class ModeExtension { Stereo = 0, IntensityStereo = 1, MsStereo = 2 };

bool operator &(ModeExtension m1, ModeExtension m2) {
    return static_cast<int>(m1) & static_cast<int>(m2);
}

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
    ModeExtension mode_extension;
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
    int samples[2][2][576] = {};
    double requantized_samples[2][2][576] = {};
    double output[2][2][32][18] = {};
};

float V1[1024] = {};
float V2[1024] = {};

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
    header.mode_extension       = static_cast<ModeExtension>(bitstream.read_bits(2));
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
    if (size == 1) {
        return *table;
    }

    int result = bitstream.read_bit();
    int bits_read = 1;
    int bits[22] = {};
    bits[0] = result;

    while (bits_read < 22 && !bitstream.eos()) {
        for (size_t i = 0; i < size; i++) {
            if (table[i].hcod == result && table[i].hlen == bits_read) {
                return table[i];
            }
        }

        result = (result << 1) | bitstream.read_bit();
        bits[bits_read] = result & 1;
        ++bits_read;
    }

    bitstream.rewind(bits_read);
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

    bitstream.rewind(bits_read);
    throw std::exception();
}

double requantize_III(int sample, double exponent) {
    const int sign = sample < 0 ? -1 : 1;
    const int magnitude = abs(sample);
    double result = sign * pow(magnitude, 4 / 3.0) * exponent;
    if (std::isnan(exponent) || std::isnan(result)) {
        std::cout << "ISNAN!" << std::endl;
    }
    return result;
}

void exponents_III(
    const Header& header, 
    const bool mixed_mode, 
    const bool has_preflag, 
    const bool scalefac_scale, 
    const int block_type, 
    const int global_gain, 
    const int *subblock_gain, 
    const int scalefac_l[22], 
    const int scalefac_s[39], 
    const ScaleFactorBand *sfb, 
    double exp[576]) {
    
    auto fill_band = [exp](double exponent, int start, int end) {
        for (int j = start; j <= end && j < 576; j++) {
            exp[j] = exponent;
        }
    };

    double scalefac_multiplier = scalefac_scale ? 1 : 0.5;
    int gain = global_gain - 210;

    if (block_type != 2) {
        for (size_t sfbi = 0; sfbi < 22; sfbi++) {
            double exponent = gain/4.0 - (scalefac_multiplier * (scalefac_l[sfbi] + has_preflag * layer_III_pretab[sfbi]));
            fill_band(pow(2, exponent), sfb[sfbi].start, sfb[sfbi].end);
        }
    } else {
        size_t sfbi = 0;
        size_t l = 0;

        if (mixed_mode) {
            // TODO: need mixed mode short scale factor band lists
            while (l < 36) {
                int exponent = gain - (scalefac_multiplier * (scalefac_l[sfbi] + has_preflag * layer_III_pretab[sfbi]));
                fill_band(pow(2, exponent), sfb[sfbi].start, sfb[sfbi].end);
                l += sfb[sfbi++].width;
            }
        }

        double gain0 = (gain - 8 * subblock_gain[0])/4.0;
        double gain1 = (gain - 8 * subblock_gain[1])/4.0;
        double gain2 = (gain - 8 * subblock_gain[2])/4.0;

        while (l < 576 && sfbi < 39) {
            double exponent0 = gain0 - (scalefac_multiplier * scalefac_s[sfbi + 0]);
            double exponent1 = gain1 - (scalefac_multiplier * scalefac_s[sfbi + 1]);
            double exponent2 = gain2 - (scalefac_multiplier * scalefac_s[sfbi + 2]);

            fill_band(pow(2, exponent0), sfb[sfbi + 0].start, sfb[sfbi + 0].end);
            l += sfb[sfbi + 0].width;
            fill_band(pow(2, exponent1), sfb[sfbi + 1].start, sfb[sfbi + 1].end);
            l += sfb[sfbi + 1].width;
            fill_band(pow(2, exponent2), sfb[sfbi + 2].start, sfb[sfbi + 2].end);
            l += sfb[sfbi + 2].width;
 
            sfbi += 3;
        }

        while (l < 576) {
            exp[l++] = 0;
        }
    }
}

void reorder_III(double samples[576], bool mixed_mode, ScaleFactorBand* sfb) {
    size_t sfbi = 0;
    size_t l = 0;
    double tmp[576] = {};

    if (mixed_mode) {
        // TODO: need mixed mode short scale factor band lists
        while (l < 36) {
            for (int i = 0; i < sfb[sfbi].width; i++) {
                tmp[l] = samples[l];
            }
            l += sfb[sfbi++].width;
        }
    }

    while (l < 576 && sfbi <= 36) {
        for (int i = 0; i < sfb[sfbi].width; i++) {
            tmp[l++] = samples[sfb[sfbi + 0].start + i];
            tmp[l++] = samples[sfb[sfbi + 1].start + i];
            tmp[l++] = samples[sfb[sfbi + 2].start + i];
        }
        sfbi += 3;
    }

    memcpy(samples, tmp, 576 * sizeof(double));
}

void alias_reduce_III(AudioDataIII& data, int ch, int gr) {
    for (size_t sb = 0; sb < 576 - 18; sb += 18) {
        for (size_t i = 0; i < 8; i++) {
            const int idx1 = sb + 17 - i;
            const int idx2 = sb + 18 + i;
            double d1 = data.requantized_samples[ch][gr][idx1];
            double d2 = data.requantized_samples[ch][gr][idx2];
            data.requantized_samples[ch][gr][idx1] = d1 * layer_III_cs[i] - d2 * layer_III_ca[i];
            data.requantized_samples[ch][gr][idx2] = d2 * layer_III_cs[i] + d1 * layer_III_ca[i];
        }
    }
}

void imdct_III(double input[18], double output[36], int block_type) {
    if (block_type == 2) {
        const int n = 12;
        double temp[36];

        for (size_t i = 0; i < n; i++) {
            temp[i+0] = 0;
            for (size_t k = 0; k < n / 2; k++) {
                temp[i+0] += input[3*k+0] * cos(M_PI / (2 * n) * (2 * i + 1 + n / 2) * (2 * k + 1));
            }
            temp[i+0] *= layer_III_window_2[i];
        }

        for (size_t i = 0; i < n; i++) {
            temp[i+12] = 0;
            for (size_t k = 0; k < n / 2; k++) {
                temp[i+12] += input[3*k+1] * cos(M_PI / (2 * n) * (2 * i + 1 + n / 2) * (2 * k + 1));
            }
            temp[i+12] *= layer_III_window_2[i];
        }

        for (size_t i = 0; i < n; i++) {
            temp[i+24] = 0;
            for (size_t k = 0; k < n / 2; k++) {
                temp[i+24] += input[3*k+2] * cos(M_PI / (2 * n) * (2 * i + 1 + n / 2) * (2 * k + 1));
            }
            temp[i+24] *= layer_III_window_2[i];
        }

        double* idmct1 = &temp[0];
        double* idmct2 = &temp[12];
        double* idmct3 = &temp[24];
        for (size_t i = 0; i < 6; i++) output[i] = 0;
        for (size_t i = 6; i < 12; i++) output[i] = idmct1[i - 6];
        for (size_t i = 12; i < 18; i++) output[i] = idmct1[i - 6] + idmct2[i - 12];
        for (size_t i = 18; i < 24; i++) output[i] = idmct2[i - 12] + idmct3[i - 18];
        for (size_t i = 24; i < 30; i++) output[i] = idmct3[i - 18];
        for (size_t i = 30; i < 36; i++) output[i] = 0;

    } else {
        const int n = 36;
        for (size_t i = 0; i < n; i++) {
            output[i] = 0;
            for (size_t k = 0; k < n / 2; k++) {
                output[i] += input[k] * cos(M_PI / (2 * n) * (2 * i + 1 + n / 2) * (2 * k + 1));
            }

            switch (block_type) {
            case 0:
                output[i] *= layer_III_window_0[i];
                break;
            case 1:
                output[i] *= layer_III_window_1[i];
                break;
            case 3:
                output[i] *= layer_III_window_3[i];
                break;
            }
        }
    }
}

struct layer_III_sideinfo {
    int main_data_begin = {};
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
    int scalefac_l[2][21] = {};
    int scalefac_s[2][39] = {};
};

void read_side_info_III(const Header& header, layer_III_sideinfo& si, BitStream& bitstream) {
    const int channels = header.mode == Mode::SingleChannel ? 1 : 2;

    si.main_data_begin = bitstream.read_bits(9);

    if (header.mode == Mode::SingleChannel) {
        si.private_bits = bitstream.read_bits(5);
    } else {
        si.private_bits = bitstream.read_bits(3);
    }

    for (size_t ch = 0; ch < channels; ch++) {
        for (size_t scfsi_band = 0; scfsi_band < 4; scfsi_band++) {
            si.scfsi[ch][scfsi_band] = bitstream.read_bits(1);
        }
    }

    for (size_t gr = 0; gr < 2; gr++) {
        for (size_t ch = 0; ch < channels; ch++) {
            si.part2_3_length[gr][ch] = bitstream.read_bits(12);
            si.big_values[gr][ch] = bitstream.read_bits(9);
            si.global_gain[gr][ch] = bitstream.read_bits(8);
            si.scalefac_compress[gr][ch] = bitstream.read_bits(4);
            si.window_switching_flag[gr][ch] = bitstream.read_bits(1);
            if (si.window_switching_flag[gr][ch]) {
                si.block_type[gr][ch] = bitstream.read_bits(2);
                si.mixed_block_flag[gr][ch] = bitstream.read_bits(1);
                for (size_t region = 0; region < 2; region++)
                    si.table_select[gr][ch][region] = bitstream.read_bits(5);
                for (size_t window = 0; window < 3; window++)
                    si.sub_block_gain[gr][ch][window] = bitstream.read_bits(3);
                si.region0_count[gr][ch] = (si.block_type[gr][ch] == 2 && !si.mixed_block_flag[gr][ch]) ? 8 : 7;
                si.region1_count[gr][ch] = (si.block_type[gr][ch] == 2 && !si.mixed_block_flag[gr][ch]) ? 36 : (20 - si.region0_count[gr][ch]);
            } else {
                for (size_t region = 0; region < 3; region++)
                    si.table_select[gr][ch][region] = bitstream.read_bits(5);
                si.region0_count[gr][ch] = bitstream.read_bits(4);
                si.region1_count[gr][ch] = bitstream.read_bits(3);
            }
            si.preflag[gr][ch] = bitstream.read_bits(1);
            si.scalefac_scale[gr][ch] = bitstream.read_bits(1);
            si.count1table_select[gr][ch] = bitstream.read_bits(1);
        }
    }
}

void read_scale_factors_III(const Header& header, layer_III_sideinfo& si, RingBitStream& reservoir, int gr, int ch) {
    if (si.window_switching_flag[gr][ch] == 1 && si.block_type[gr][ch] == 2) {
        if (si.mixed_block_flag[gr][ch]) {
            for (size_t sfb = 0; sfb < 8; sfb++) {
                const int bits = layer_III_scalefac_compress_slen1[si.scalefac_compress[gr][ch]];
                si.scalefac_l[ch][sfb] = reservoir.read_bits(bits);
            }
            for (size_t sfb = 3; sfb < 12; sfb++) {
                const int bits = sfb <= 5 ? layer_III_scalefac_compress_slen1[si.scalefac_compress[gr][ch]] :
                                            layer_III_scalefac_compress_slen2[si.scalefac_compress[gr][ch]];
                si.scalefac_s[ch][3 * sfb + 0] = reservoir.read_bits(bits);
                si.scalefac_s[ch][3 * sfb + 1] = reservoir.read_bits(bits);
                si.scalefac_s[ch][3 * sfb + 2] = reservoir.read_bits(bits); 
            }
        } else {
            for (size_t sfb = 0; sfb < 12; sfb++) {
                const int bits = sfb <= 5 ? layer_III_scalefac_compress_slen1[si.scalefac_compress[gr][ch]] :
                                            layer_III_scalefac_compress_slen2[si.scalefac_compress[gr][ch]];
                si.scalefac_s[ch][3 * sfb + 0] = reservoir.read_bits(bits);
                si.scalefac_s[ch][3 * sfb + 1] = reservoir.read_bits(bits);
                si.scalefac_s[ch][3 * sfb + 2] = reservoir.read_bits(bits); 
            }
        }
    } else {
        auto scalefactor_bits = [&gr, &ch, &si](int sfb) -> int {
            return sfb <= 10 ? layer_III_scalefac_compress_slen1[si.scalefac_compress[gr][ch]]
                             : layer_III_scalefac_compress_slen2[si.scalefac_compress[gr][ch]];
        };

        if ((si.scfsi[ch][0] == 0) || (gr == 0)) {
            for (size_t sfb = 0; sfb < 6; sfb++) {
                const int bits = scalefactor_bits(sfb);
                si.scalefac_l[ch][sfb] = reservoir.read_bits(bits);
            }
        }
        if ((si.scfsi[ch][1] == 0) || (gr == 0)) {
            for (size_t sfb = 6; sfb < 11; sfb++) {
                const int bits = scalefactor_bits(sfb);
                si.scalefac_l[ch][sfb] = reservoir.read_bits(bits);
            }
        }
        if ((si.scfsi[ch][2] == 0) || (gr == 0)) {
            for (size_t sfb = 11; sfb < 16; sfb++) {
                const int bits = scalefactor_bits(sfb);
                si.scalefac_l[ch][sfb] = reservoir.read_bits(bits);
            }
        }
        if ((si.scfsi[ch][3] == 0) || (gr == 0)) {
            for (size_t sfb = 16; sfb < 21; sfb++) {
                const int bits = scalefactor_bits(sfb);
                si.scalefac_l[ch][sfb] = reservoir.read_bits(bits);
            }
        }
    }
}

int read_huffman_data_III(const Header& header,
                          const layer_III_sideinfo& si,
                          RingBitStream& reservoir,
                          int gr, int ch,
                          AudioDataIII& result,
                          int granule_bits_read) {

    double exponents[576] = {};
    exponents_III(
        header,
        si.mixed_block_flag[gr][ch],
        si.preflag[gr][ch],
        si.scalefac_scale[gr][ch],
        si.block_type[gr][ch],
        si.global_gain[gr][ch],
        si.sub_block_gain[gr][ch],
        si.scalefac_l[ch],
        si.scalefac_s[ch],
        si.block_type[gr][ch] == 2 ? ScaleFactorBandsShort[header.sampling_frequency].data() :
                                     ScaleFactorBandsLong[header.sampling_frequency].data(),
        exponents);

    int count = 0;
    int scaleFactorBandIndex1 = si.region0_count[gr][ch] + 1;
    int scaleFactorBandIndex2 = scaleFactorBandIndex1 + si.region1_count[gr][ch] + 1;
    if (scaleFactorBandIndex2 >= ScaleFactorBandsLong[header.sampling_frequency].size()) {
        scaleFactorBandIndex2 = ScaleFactorBandsLong[header.sampling_frequency].size() - 1;
    }

    int region1_start = ScaleFactorBandsLong[header.sampling_frequency][scaleFactorBandIndex1].start;
    int region2_start = ScaleFactorBandsLong[header.sampling_frequency][scaleFactorBandIndex2].start;
    if (si.window_switching_flag[gr][ch] && si.block_type[gr][ch] == 2) {
        region1_start = 36;
        region2_start = 576;
    }

    const HuffmanTable2* table = nullptr;
    size_t i = 0;
    for (; i < si.big_values[gr][ch] * 2; i += 2) {
        if (i < region1_start) {
            table = &HuffmanTables2[si.table_select[gr][ch][0]];
        } else if (i < region2_start) {
            table = &HuffmanTables2[si.table_select[gr][ch][1]];
        } else {
            table = &HuffmanTables2[si.table_select[gr][ch][2]];
        }

        if (table == nullptr || table->table == nullptr) {
            std::cout << "UNKNOWN TABLE: " << si.table_select[gr][ch][0] << std::endl;
            throw std::exception();
        }

        const HuffmanEntry2& huff = huffman_decode(reservoir, table->table, table->table_size);
        int x = huff.x;
        int y = huff.y;
        granule_bits_read += huff.hlen;
        
        if (x == 15 && table->linbits > 0) {
            x += reservoir.read_bits(table->linbits);
            granule_bits_read += table->linbits;
        }
        if (x != 0) {
            if (reservoir.read_bit()) x = -x;
            granule_bits_read++;
        }

        if (y == 15 && table->linbits > 0) {
            y += reservoir.read_bits(table->linbits);
            granule_bits_read += table->linbits;
        }
        if (y != 0) {
            if (reservoir.read_bit()) y = -y;
            granule_bits_read++;
        }

        result.samples[ch][gr][count+0] = x;
        result.samples[ch][gr][count+1] = y;

        result.requantized_samples[ch][gr][count + 0] = requantize_III(x, exponents[count + 0]);
        result.requantized_samples[ch][gr][count + 1] = requantize_III(y, exponents[count + 1]);

        count += 2;
    }

    const HuffmanEntry4* count1table = si.count1table_select[gr][ch] ? huffman_table_B : huffman_table_A;

    // count1 is not known. We have to read huffman encoded values
    // until we've exhausted the granule's bits. We know the size of
    // the granule from part2_3_length, which is the number of bits
    // used for scaleactors and huffman data (in the granule).
    while (granule_bits_read < si.part2_3_length[gr][ch] && count < 576) {
        const HuffmanEntry4 *entry;
        try {
            entry = &huffman_decode(reservoir, count1table, 16);
        } catch (std::exception e) {
            break;
        }
        int v = entry->v;
        if (v != 0) {
            if (reservoir.read_bit()) v = -v;
            granule_bits_read++;
        }
        int w = entry->w;
        if (w != 0) {
            if (reservoir.read_bit()) w = -w;
            granule_bits_read++;
        }
        int x = entry->x;
        if (x != 0) {
            if (reservoir.read_bit()) x = -x;
            granule_bits_read++;
        }
        int y = entry->y;
        if (y != 0) {
            if (reservoir.read_bit()) y = -y;
            granule_bits_read++;
        }

        result.samples[ch][gr][count+0] = v;
        result.samples[ch][gr][count+1] = w;
        result.samples[ch][gr][count+2] = x;
        result.samples[ch][gr][count+3] = y;

        result.requantized_samples[ch][gr][count + 0] = requantize_III(v, exponents[count + 0]);
        result.requantized_samples[ch][gr][count + 1] = requantize_III(w, exponents[count + 1]);
        result.requantized_samples[ch][gr][count + 2] = requantize_III(x, exponents[count + 2]);
        result.requantized_samples[ch][gr][count + 3] = requantize_III(y, exponents[count + 3]);

        count += 4;

        granule_bits_read += entry->hlen;
    }

    if (granule_bits_read > si.part2_3_length[gr][ch]) {
        reservoir.rewind(granule_bits_read - si.part2_3_length[gr][ch]);
        count--;
    }

    if (granule_bits_read < si.part2_3_length[gr][ch]) {
        for (size_t i = granule_bits_read; i < si.part2_3_length[gr][ch] && !reservoir.eos(); i++) {
            reservoir.read_bit();
        }
    }

    return count;
}

int get_first_zero_band(double* samples, ScaleFactorBand* bands, size_t size) {
    int last_nonempty_band = 0;

    for (size_t i = 0; i < size; i++) {
        bool is_empty = true;
        for (size_t l = bands[i].start; l < bands[i].end; l++) {
            if (samples[l] != 0) {
                is_empty = false;
                break;
            }
        }
        if (!is_empty) {
            last_nonempty_band = i;
        }
    }

    return last_nonempty_band;
}

void stereo_III(const Header& header, AudioDataIII& data, int gr, const layer_III_sideinfo& si) {
    const double SQRT_2 = 1.4142135623730950488016887242097;

    size_t sfbi_ms_start = 0;
    size_t sfbi_ms_end = 0;
    size_t sfbi_intensity_start = 0;
    size_t sfbi_intensity_end = 0;
    ScaleFactorBand* sfbs = nullptr;
    size_t sfbs_length = 0;

    if (si.block_type[gr][0] == 2) {
        sfbs = ScaleFactorBandsShort[header.sampling_frequency].data();
        sfbs_length = ScaleFactorBandsShort[header.sampling_frequency].size();
    } else {
        sfbs = ScaleFactorBandsLong[header.sampling_frequency].data();
        sfbs_length = ScaleFactorBandsLong[header.sampling_frequency].size();
    }

    if (header.mode_extension & ModeExtension::IntensityStereo) {
        sfbi_intensity_start = get_first_zero_band(data.requantized_samples[1][gr], sfbs, sfbs_length);
        sfbi_intensity_end = sfbs_length;
    }

    if (header.mode_extension & ModeExtension::MsStereo) {
        sfbi_ms_start = 0;
        sfbi_ms_end = (header.mode_extension & ModeExtension::IntensityStereo) ? sfbi_intensity_start : sfbs_length;
    }

    for (size_t sfbi = sfbi_ms_start; sfbi < sfbi_ms_end; sfbi++) {
        for (size_t i = sfbs[sfbi].start; i <= sfbs[sfbi].end; i++) {
            const double m = data.requantized_samples[0][gr][i];
            const double s = data.requantized_samples[1][gr][i];
            data.requantized_samples[0][gr][i] = (m + s) / SQRT_2;
            data.requantized_samples[1][gr][i] = (m - s) / SQRT_2;
        }
    }

    for (size_t sfbi = sfbi_intensity_start; sfbi < sfbi_intensity_end; sfbi++) {
        auto is_pos = si.scalefac_l[1][sfbi];
        if (is_pos >= 7) continue;
        double is_ratio = tan(is_pos * M_PI / 12);
        for (size_t i = sfbs[sfbi].start; i <= sfbs[sfbi].end; i++) {
            double sample_left = data.requantized_samples[0][gr][i];
            double coeff_l = is_ratio / (1 + is_ratio);
            double coeff_r = 1 / (1 + is_ratio);
            data.requantized_samples[0][gr][i] = sample_left * coeff_l;
            data.requantized_samples[1][gr][i] = sample_left * coeff_r;
        }
    }
}

AudioDataIII read_audio_data_III(const Header& header, BitStream& bitstream, RingBitStream& reservoir, double lastValues[2][32][18]) {
    layer_III_sideinfo si;
    const int channels = header.mode == Mode::SingleChannel ? 1 : 2;

    read_side_info_III(header, si, bitstream);

    std::cout << "Gr0 bt " << si.block_type[0][0] << ", Gr1 bt " << si.block_type[1][0] << std::endl;

    // Copy all the remaining data in this frame into the bit reservoir.
    // Append it to the left-over data of the last frame in order to build the
    // complete current frame.
    {
        _ASSERT(bitstream.get_current_bit() == 0);
        reservoir.seek_to_end();

        for (size_t i = (channels == 2 ? 32 : 17) + 4; i < header.frame_size; i++) {
            unsigned char data = bitstream.read_bits<unsigned char>(8);
            reservoir.write(&data, 1);
        }

        if (si.main_data_begin > 0) {
            reservoir.seek_relative(-si.main_data_begin);
        }
    }

    AudioDataIII result;

    std::fstream f;
    f.open("out.me", std::fstream::out | std::fstream::app);

    for (size_t gr = 0; gr < 2; gr++) {
        for (size_t ch = 0; ch < channels; ch++) {
            const int current_position = reservoir.position();
            read_scale_factors_III(header, si, reservoir, gr, ch);
            read_huffman_data_III(header, si, reservoir, gr, ch, result, reservoir.position() - current_position);

            f << "gr " << gr << " ch " << ch << std::endl;
            for (size_t i = 0; i < 32; i++) {
                for (size_t j = 0; j < 18; j++) {
                    //f << result.requantized_samples[ch][gr][i * 18 + j] << ", ";
                    f << std::round(result.requantized_samples[ch][gr][i * 18 + j] * 100000000)/100000000.0 << ", ";
                }
                f << std::endl;
            }
            f << std::endl;

            if (si.block_type[gr][ch] == 2) {
                reorder_III(result.requantized_samples[ch][gr], si.mixed_block_flag[gr][ch], ScaleFactorBandsShort[header.sampling_frequency].data());
            } else {
                alias_reduce_III(result, ch, gr);
            }
        }
        stereo_III(header, result, gr, si);
    }

    for (size_t gr = 0; gr < 2; gr++) {
        for (size_t ch = 0; ch < channels; ch++) {
            for (size_t i = 0; i < 576; i += 18) {
                double output[36];
                imdct_III(&result.requantized_samples[ch][gr][i], output, si.block_type[gr][ch]);

                const int sb = i / 18;
                for (size_t s = 0; s < 18; s++) {
                    // overlap add
                    result.output[gr][ch][sb][s] = output[s] + lastValues[ch][sb][s];
                    lastValues[ch][sb][s] = output[s + 18];

                    // frequency inversion
                    if (sb % 2 == 1 && s % 2 == 1) result.output[gr][ch][sb][s] *= -1;
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
    const int bound = header.mode == Mode::IntensityStereo ? table_bounds[static_cast<int>(header.mode_extension)] : 32;

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
void synthesis(float *V, float samples[32], float result[32]) {
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
        std::ifstream infile(R"(..\res\testcase.mp3)", std::ios::binary);
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

    float in_samples[2][32];
    float out_samples[2][32];
    int channels = 0;

    auto write_samples = [&channels, &out_samples, &outfile, &sample_count]() {
        for (size_t j = 0; j < 32; j++) {
            const short sample_l = clamp(32767 * out_samples[0][j]);
            outfile.write(reinterpret_cast<const char*>(&sample_l), 2);
            if (channels == 2) {
                const short sample_r = clamp(32767 * out_samples[1][j]);
                outfile.write(reinterpret_cast<const char*>(&sample_r), 2);
            }
        }
        sample_count += 32;
    };

    RingBitStream reservoir { 65536 };
    double lastValues[2][32][18] = {};
    int frameCount = 0;

    while (!bitstream.eof()) {
        if (!synchronize(bitstream)) {
            std::cout << "no sync found!" << std::endl;
            break;
        }
        std::cout << "Found frame @ " << bitstream.get_current_byte() << std::endl;

        const size_t frame_position = bitstream.get_current_byte() - 1;

        Header header = read_header(bitstream);
        std::cout << "Layer " << header.layer << ", " << header.bitrate << ", " << header.padding_bit << ", " << header.protection_bit << ", " << header.frame_size << " (" << frameCount << ")" << std::endl;

        if (header.layer != 3) {
            std::cout << "ERRROR: LOST SYNC!" << std::endl;
            continue;
        }

        ++frameCount;
        if (header.layer == 1) {
            AudioDataI audioData = read_audio_data_I(header, bitstream);

            for (size_t i = 0; i < 12; i++) {
                for (size_t j = 0; j < 32; j++) {
                    in_samples[0][j] = audioData.requantized_samples[0][j][i];
                }
                synthesis(V1, in_samples[0], out_samples[0]);
                write_samples();
            }
        } else if (header.layer == 2) {
            AudioDataII audioData = read_audio_data_II(header, bitstream);

            for (size_t i = 0; i < 36; i++) {
                for (size_t j = 0; j < 32; j++) {
                    in_samples[0][j] = audioData.requantized_samples[0][j][i];
                }
                synthesis(V1, in_samples[0], out_samples[0]);
                write_samples();
            }
        } else if (header.layer == 3) {
            AudioDataIII audioData = read_audio_data_III(header, bitstream, reservoir, lastValues);
            for (size_t gr = 0; gr < 2; gr++) {
                for (size_t i = 0; i < 18; i++) {
                    for (size_t j = 0; j < 32; j++) {
                        in_samples[0][j] = 1 * audioData.output[gr][0][j][i];
                        in_samples[1][j] = 1 * audioData.output[gr][1][j][i];
                    }
                    synthesis(V1, in_samples[0], out_samples[0]);
                    synthesis(V2, in_samples[1], out_samples[1]);
                    write_samples();
                }
            }
        }
        else {
            bitstream.go_to_byte(bitstream.get_current_byte() - 2);
        }
    }

    write_wav_header(outfile, 44100, channels, sample_count);
}
