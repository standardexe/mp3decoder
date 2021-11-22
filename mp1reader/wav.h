#pragma once

#include <fstream>

typedef struct WAV_HEADER {
    unsigned char RIFF[4] = { 'R', 'I', 'F', 'F' };
    unsigned int ChunkSize;
    unsigned char WAVE[4] = { 'W', 'A', 'V', 'E' };
    unsigned char fmt[4] = { 'f', 'm', 't', ' ' };
    unsigned int Subchunk1Size = 16;
    unsigned short AudioFormat = 1;
    unsigned short NumOfChan = 1;
    unsigned int SamplesPerSec = 16000;
    unsigned int bytesPerSec = 16000 * 2;
    unsigned short blockAlign = 2;
    unsigned short bitsPerSample = 16;
    unsigned char Subchunk2ID[4] = { 'd', 'a', 't', 'a' };
    unsigned int Subchunk2Size;
} wav_hdr;

void write_wav_header(std::ofstream& file, int samplerate, int channels, size_t sample_count) {
    wav_hdr hdr;
    hdr.SamplesPerSec = samplerate;
    hdr.NumOfChan = channels;
    hdr.bitsPerSample = 16;
    hdr.bytesPerSec = hdr.SamplesPerSec * hdr.NumOfChan * hdr.bitsPerSample / 8;
    hdr.blockAlign = hdr.NumOfChan * (hdr.bitsPerSample / 8);
    hdr.Subchunk2Size = sample_count * hdr.blockAlign;
    hdr.ChunkSize = hdr.Subchunk2Size + 44;

    file.seekp(0, std::ios::_Seekbeg);
    file.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
}