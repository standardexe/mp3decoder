#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "bitbuf_view.h"
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
    int scalefactors[2][32] = {};
    int samples[2][32][12] = {};
    float requantized_samples[2][32][12] = {};
};

struct AudioDataII {
    int allocations[2][32] = {};
    int scfsi[2][32] = {};
    float scalefactors[2][32][3] = {};
    int samples[2][32][36] = {};
    float requantized_samples[2][32][36] = {};
};

float V[1024] = {};

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
    header.layer                = table_layer[bitstream.read_bits<>(2)];
    header.protection_bit       = bitstream.read_bit();
    header.bitrate              = table_bitrate_per_layer[header.layer - 1][bitstream.read_bits<>(4)];
    header.sampling_frequency   = table_samplerate[bitstream.read_bits<>(2)];
    header.padding_bit          = bitstream.read_bit();
    header.private_bit          = bitstream.read_bit();
    header.mode                 = static_cast<Mode>(bitstream.read_bits<>(2));
    header.mode_extension       = bitstream.read_bits<>(2);
    header.copyright_bit        = bitstream.read_bit();
    header.original_bit         = bitstream.read_bit();
    header.emphasis             = static_cast<Emphasis>(bitstream.read_bits<>(2));
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

size_t left_most_bit_index(int value) {
    size_t counter = 0;
    while (value > 0) {
        value >>= 1;
        ++counter;
    }
    return counter;
}

AudioDataII read_audio_data_II(const Header& header, BitStream& bitstream) {
    AudioDataII data;

    auto quantization_info = layer_II_get_quantization_table(header);

    const int channels = header.mode == Mode::SingleChannel ? 1 : 2;
    const int sblimit = quantization_info.sblimit;
    const int bound = sblimit; // TODO: for Joint Stereo this is not true

    for (size_t sb = 0; sb < bound; sb++) {
        for (size_t ch = 0; ch < channels; ch++) {
            const int nbal = (*quantization_info.nbal)[sb];
            const int al = bitstream.read_bits<>(nbal);
            assert(al >= 0 && al <= 15);
            data.allocations[ch][sb] = (*quantization_info.table)[sb][al];
            assert(data.allocations[ch][sb] > -1);
        }
    }

    for (size_t sb = bound; sb < sblimit; sb++) {
        const int al = bitstream.read_bits<>((*quantization_info.nbal)[sb]);
        data.allocations[0][sb]  = bitstream.read_bits<>((*quantization_info.nbal)[sb]);
        data.allocations[1][sb] = data.allocations[0][sb];
        assert(data.allocations[0][sb] > -1);
    }

    for (size_t sb = 0; sb < sblimit; sb++) {
        for (size_t ch = 0; ch < channels; ch++) {
            if (data.allocations[ch][sb] != 0) {
                data.scfsi[ch][sb] = bitstream.read_bits<>(2);
            }
        }
    }

    for (size_t sb = 0; sb < sblimit; sb++) {
        for (size_t ch = 0; ch < channels; ch++) {
            if (data.allocations[ch][sb] != 0) {
                switch (data.scfsi[ch][sb]) {
                case 0:
                    data.scalefactors[ch][sb][0] = table_scale_factors[bitstream.read_bits<>(6)];
                    data.scalefactors[ch][sb][1] = table_scale_factors[bitstream.read_bits<>(6)];
                    data.scalefactors[ch][sb][2] = table_scale_factors[bitstream.read_bits<>(6)];
                    break;
                case 2:
                    data.scalefactors[ch][sb][0] = table_scale_factors[bitstream.read_bits<>(6)];
                    data.scalefactors[ch][sb][1] = data.scalefactors[ch][sb][0];
                    data.scalefactors[ch][sb][2] = data.scalefactors[ch][sb][0];
                    break;
                case 1:
                    data.scalefactors[ch][sb][0] = table_scale_factors[bitstream.read_bits<>(6)];
                    data.scalefactors[ch][sb][1] = data.scalefactors[ch][sb][0];
                    data.scalefactors[ch][sb][2] = table_scale_factors[bitstream.read_bits<>(6)];
                    break;
                case 3:
                    data.scalefactors[ch][sb][0] = table_scale_factors[bitstream.read_bits<>(6)];
                    data.scalefactors[ch][sb][1] = table_scale_factors[bitstream.read_bits<>(6)];
                    data.scalefactors[ch][sb][2] = data.scalefactors[ch][sb][1];
                    break;
                }
            }
        }
    }

    auto quantization_class_index = [](int value) -> int {
        for (size_t i = 0; i < 17; i++) {
            if (layer_II_quantization_class_num_steps[i] == value) {
                return i;
            }
        }
        throw std::exception();
    };

    auto grouping = [](int quant_index) -> bool {
        return layer_II_quantization_class_group[quant_index];
    };

    for (size_t gr = 0; gr < 12; gr++) {
        for (size_t sb = 0; sb < bound; sb++) {
            for (size_t ch = 0; ch < channels; ch++) {
                if (data.allocations[ch][sb] != 0) {
                    const int quant_index = quantization_class_index(data.allocations[ch][sb]);
                    const int code_width = layer_II_quantization_class_bits_per_cw[quant_index];
                    if (grouping(quant_index)) {
                        assert(code_width >= 5 && code_width <= 10);
                        const int nlevels = data.allocations[ch][sb];
                        const size_t lmbi = left_most_bit_index(nlevels);
                        int samplecode = bitstream.read_bits<>(code_width);
                        for (size_t i = 0; i < 3; i++) {
                            data.samples[ch][sb][3 * gr + i] = samplecode % nlevels;
                            samplecode = samplecode / nlevels;

                            const float scalefactor = data.scalefactors[ch][sb][gr / 4];
                            const int sample = data.samples[ch][sb][3 * gr + i];
                            const float requantized = requantize_II(sample, lmbi, quant_index);

                            data.requantized_samples[ch][sb][3 * gr + i] = scalefactor * requantized;
                        }
                    } else {
                        assert(code_width >= 3 && code_width <= 16);
                        for (size_t s = 0; s < 3; s++) {
                            data.samples[ch][sb][3 * gr + s] = bitstream.read_bits<>(code_width);
                            data.requantized_samples[ch][sb][3 * gr + s] = data.scalefactors[ch][sb][gr / 4] * requantize_II(data.samples[ch][sb][3 * gr + s], code_width, quant_index);
                        }
                    }
                }
            }
        }
        for (size_t sb = bound; sb < sblimit; sb++) {
            if (data.allocations[0][sb] != 0) {
                const int quant_index = quantization_class_index(data.allocations[0][sb]);
                const int code_width = layer_II_quantization_class_bits_per_cw[quant_index];
                if (grouping(quant_index)) {
                    assert(code_width >= 5 && code_width <= 10);
                    const int nlevels = data.allocations[0][sb];
                    const size_t lmbi = left_most_bit_index(nlevels);
                    int samplecode = bitstream.read_bits<>(code_width);
                    for (size_t i = 0; i < 3; i++) {
                        data.samples[0][sb][3 * gr + i] = samplecode % nlevels;
                        samplecode /= nlevels;
                        data.requantized_samples[0][sb][3 * gr + i] = data.scalefactors[0][sb][i] * requantize_II(data.samples[0][sb][3 * gr + i], lmbi, quant_index);
                    }
                } else {
                    assert(code_width >= 3 && code_width <= 16);
                    for (size_t s = 0; s < 3; s++) {
                        data.samples[0][sb][3 * gr + s] = bitstream.read_bits<>(code_width);
                        data.requantized_samples[0][sb][3 * gr + s] = data.scalefactors[0][sb][s] * requantize_II(data.samples[0][sb][3 * gr + s], code_width, quant_index);
                    }
                }
            }
        }
    }

    return data;
}

AudioDataI read_audio_data_I(Header& header, BitStream& bitstream) {
    AudioDataI data;

    const int channels = header.mode == Mode::SingleChannel ? 1 : 2;
    const int bound = header.mode == Mode::IntensityStereo ? table_bounds[header.mode_extension] : 32;

    for (size_t sb = 0; sb < bound; sb++) {
        for (size_t ch = 0; ch < channels; ch++) {
            data.allocations[ch][sb] = table_allocation[bitstream.read_bits<>(4)];
        }
    }

    for (size_t sb = bound; sb < 32; sb++) {
        data.allocations[0][sb] = bitstream.read_bits<>(4);
        data.allocations[1][sb] = data.allocations[0][sb];
    }

    for (size_t sb = 0; sb < 32; sb++) {
        for (size_t ch = 0; ch < channels; ch++) {
            if (data.allocations[ch][sb] != 0) {
                data.scalefactors[ch][sb] = bitstream.read_bits<>(6);
            }
        }
    }

    for (size_t s = 0; s < 12; s++) {
        for (size_t sb = 0; sb < bound; sb++) {
            for (size_t ch = 0; ch < channels; ch++) {
                if (data.allocations[ch][sb] != 0) {
                    const int num_bits = data.allocations[ch][sb];
                    const float scale_factor = table_scale_factors[data.scalefactors[ch][sb]];
                    data.samples[ch][sb][s] = bitstream.read_bits<>(num_bits);
                    data.requantized_samples[ch][sb][s] = scale_factor * requantize_I(data.samples[ch][sb][s], num_bits);
                }
            }
        }

        for (size_t sb = bound; sb < 32; sb++) {
            if (data.allocations[0][sb] != 0) {
                const int num_bits = data.allocations[0][sb];
                const float scale_factor = table_scale_factors[data.scalefactors[0][sb]];
                data.samples[0][sb][s] = bitstream.read_bits<>(num_bits);
                data.requantized_samples[0][sb][s] = scale_factor * requantize_I(data.samples[0][sb][s], num_bits);
            }
        }
    }

    return data;
}

void synthesis(float samples[32], float result[32]) {
    // 1. Shifting
    for (size_t i = 1023; i >= 64; i--) {
        V[i] = V[i - 64];
    }

    // 2. Matrixing
    for (size_t i = 0; i < 64; i++) {
        V[i] = 0;
        for (size_t k = 0; k < 32; k++) {
            const float N = cos((16 + i) * (2 * k + 1) * M_PI / 64.0f);
            V[i] += N * samples[k];
        }
    }

    // 3. Build values vector U
    float U[512];
    for (size_t i = 0; i < 8; i++) {
        for (size_t j = 0; j < 32; j++) {
            U[i * 64 + j] = V[i * 128 + j];
            U[i * 64 + 32 + j] = V[i * 128 + 96 + j];
        }
    }

    // 4. Window
    float W[512];
    for (size_t i = 0; i < 512; i++) {
        W[i] = U[i] * table_window[i];
    }

    // 5. Calculate 32 samples
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
        std::ifstream infile(R"(C:\Users\arne\Downloads\wav2mp1\chirp_l2.mp2)", std::ios::binary);
        infile.seekg(0, std::ios::end);
        content_length = infile.tellg();
        infile.seekg(0, std::ios::beg);

        content = new char[content_length];
        infile.read(content, content_length);
    }

    auto bitstream = BitStream(content, content_length);

    std::ofstream outfile;
    //outfile.open(R"(C:\Users\arne\Downloads\wav2mp1\out.csv)");
    outfile.open(R"(C:\Users\arne\Downloads\wav2mp1\out.wav)", std::ios::binary);
    outfile.write(new char[44], 44);
    size_t sample_count = 0;

    while (!bitstream.eof()) {
        if (!synchronize(bitstream)) {
            std::cout << "no sync found!" << std::endl;
            break;
        }

        Header header = read_header(bitstream);
        AudioDataII audioData = read_audio_data_II(header, bitstream);
        //AudioDataI audioData = read_audio_data_I(header, bitstream);

        float in_samples[32];
        float out_samples[32];
        for (size_t i = 0; i < 36; i++) {
            for (size_t j = 0; j < 32; j++) {
                in_samples[j] = audioData.requantized_samples[0][j][i];
                //outfile << in_samples[j] << ",";
            }
            //outfile << std::endl;
            synthesis(in_samples, out_samples);

            for (size_t j = 0; j < 32; j++) {
                const short sample = clamp(32767 * out_samples[j]);
                outfile.write(reinterpret_cast<const char*>(&sample), 2);
                sample_count++;
            }
        }
    }

    write_wav_header(outfile, 44100, sample_count);
}
