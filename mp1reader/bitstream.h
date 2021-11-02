#pragma once

class RingBitStream {
public:
    RingBitStream(size_t size) : size(size), write_pos(0), read_pos(0), current_bit(0) {
        data = new unsigned char[size];
    }

    void write(void* new_data, size_t new_data_size) {
        for (size_t i = 0; i < new_data_size; i++) {
            data[write_pos] = static_cast<unsigned char*>(new_data)[i];
            write_pos = (write_pos + 1) % size;
        }
    }

    bool read_bit() {
        normalize_counters();
        if (read_pos == write_pos) {
            throw std::exception(); // EOF
        }
        return (data[read_pos] >> (7 - current_bit++)) & 1;
    }

    template<typename T = int>
    T read_bits(size_t count) {
        T result = 0;
        while (count-- > 0) {
            result = result << 1 | read_bit();
        }
        return result;
    }

    void seek_to_end() {
        read_pos = write_pos;
    }

    void rewind(int bits) {
        if (current_bit > bits) {
            current_bit = 8 - bits;
        } else {
            int old_current_bit = current_bit;
            seek_relative(-1);
            current_bit = 8 - (bits - old_current_bit);
        }
    }

    void seek_relative(int bytes) {
        current_bit = 0;
        int new_pos = read_pos + bytes;
        if (new_pos < 0) {
            read_pos = size + new_pos;
        } else {
            read_pos = new_pos % size;
        }
    }

    size_t position() {
        if (read_pos > write_pos) {
            return (read_pos - write_pos) * 8 + current_bit;
        } else {
            return (size - write_pos + read_pos) * 8 + current_bit;
        }
    }

private:
    unsigned char* data;
    size_t size;
    size_t write_pos;
    size_t read_pos;
    size_t current_bit;

    void inline normalize_counters() {
        read_pos = (read_pos + current_bit / 8) % size;
        current_bit %= 8;
    }
};

class BitStream {
public:
    BitStream(void* data, size_t size) : data(static_cast<unsigned char*>(data)), size(size) {}

    bool read_bit() {
        normalize_counters();
        return (data[current_byte] >> (7 - current_bit++)) & 1;
    }

    bool eof() const { return current_byte >= size; }

    template<typename T = int>
    T read_bits(size_t count) {
        T result = 0;
        while (count-- > 0) {
            result = result << 1 | read_bit();
        }
        return result;
    }

    void go_to_next_byte() {
        if (current_bit > 0) {
            current_byte++;
            current_bit = 0;
        }
    }

    void go_to_byte(size_t pos) {
        current_bit = 0;
        current_byte = pos;
    }

    size_t get_current_bit() const { return current_bit % 8; }
    size_t get_current_byte() const { return current_byte; }
    size_t position() const { return current_byte * 8 + current_bit; }

private:
    unsigned char* data { nullptr };
    size_t size { 0 };

    size_t current_byte { 0 };
    size_t current_bit{ 0 };

    void inline normalize_counters() {
        current_byte += current_bit / 8;
        current_bit %= 8;
    }

};
