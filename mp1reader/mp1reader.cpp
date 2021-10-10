#include <iostream>
#include <fstream>
#include <vector>
#include <functional>

#include "bitstream.h"
#include "tables.h"
#include "wav.h"

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
        std::ifstream infile(R"(..\res\chirp_l1.mp2)", std::ios::binary);
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

    while (!bitstream.eof()) {
        if (!synchronize(bitstream)) {
            std::cout << "no sync found!" << std::endl;
            break;
        }

        Header header = read_header(bitstream);
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
            std::cout << "Layer III not supported" << std::endl;
            break;
        }
    }

    write_wav_header(outfile, 44100, sample_count);
}
