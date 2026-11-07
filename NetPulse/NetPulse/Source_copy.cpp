#include <iostream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <cstdint>

#ifdef _WIN32
#include <windows.h>
#endif

std::uint64_t ToUint64(const FILETIME& ft) {
	ULARGE_INTEGER value{};

	value.LowPart = ft.dwLowDateTime;
	value.HighPart = ft.dwHighDateTime;

	return value.QuadPart;
}

struct CpuTimes {
	std::uint64_t idle;
	std::uint64_t kernel;
	std::uint64_t user;
};