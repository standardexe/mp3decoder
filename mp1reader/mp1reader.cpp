#include <iostream>
#include <fstream>
#include <vector>

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
    int scalefactors[2][32] = {};
    int samples[2][32][12] = {};
    float requantized_samples[2][32][12] = {};
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

float requantize(int raw_sample, int num_bits) {
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
                    data.requantized_samples[ch][sb][s] = scale_factor * requantize(data.samples[ch][sb][s], num_bits);
                }
            }
        }

        for (size_t sb = bound; sb < 32; sb++) {
            if (data.allocations[0][sb] != 0) {
                const int num_bits = data.allocations[0][sb];
                const float scale_factor = table_scale_factors[data.scalefactors[0][sb]];
                data.samples[0][sb][s] = bitstream.read_bits<>(num_bits);
                data.requantized_samples[0][sb][s] = scale_factor * requantize(data.samples[0][sb][s], num_bits);
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
        std::ifstream infile(R"(C:\Users\arne\Downloads\wav2mp1\dolby.mp2)", std::ios::binary);
        infile.seekg(0, std::ios::end);
        content_length = infile.tellg();
        infile.seekg(0, std::ios::beg);

        content = new char[content_length];
        infile.read(content, content_length);
    }

    auto bitstream = BitStream(content, content_length);

    std::ofstream outfile;
    outfile.open(R"(C:\Users\arne\Downloads\wav2mp1\out.wav)", std::ios::binary);
    outfile.write(new char[44], 44);
    size_t sample_count = 0;

    while (!bitstream.eof()) {
        if (!synchronize(bitstream)) {
            std::cout << "no sync found!" << std::endl;
            break;
        }

        Header header = read_header(bitstream);
        AudioDataI audioData = read_audio_data_I(header, bitstream);

        float in_samples[32];
        float out_samples[32];
        for (size_t i = 0; i < 12; i++) {
            for (size_t j = 0; j < 32; j++) {
                in_samples[j] = audioData.requantized_samples[0][j][i];
            }
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
