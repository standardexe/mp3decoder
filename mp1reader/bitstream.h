#pragma once

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

	size_t get_current_bit() const { return current_bit; }
	size_t get_current_byte() const { return current_byte; }

private:
	unsigned char* data { nullptr };
	size_t size { 0 };

	size_t current_byte { 0 };
	size_t current_bit{ 0 };

	void normalize_counters() {
		current_byte += current_bit / 8;
		current_bit %= 8;
	}

};
