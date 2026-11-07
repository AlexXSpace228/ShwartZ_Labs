#include <iostream>
#include <thread>
#include <chrono>
#include <iomanip>
#include <cstdint>

#ifdef _WIN32
#include <windows.h>
#endif

std::uint64_t toUint64(const FILETIME& ft) {
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

bool getCpuTimes(CpuTimes& times) {

	FILETIME idle{}, kernel{}, user{};

	if (!GetSystemTimes(&idle, &kernel, &user)) {
		return false;
	}

	times.idle = toUint64(idle);
	times.kernel = toUint64(kernel);
	times.user = toUint64(user);

	return true;
}

int main() {

	CpuTimes previous{};

	if (!getCpuTimes(previous)) {
		std::cerr << "Failed to get CPU times\n";
		return 1;
	}

	std::cout << std::fixed << std::setprecision(1);

	while (true) {
		std::this_thread::sleep_for(std::chrono::seconds(1));

		CpuTimes current{};

		if (!getCpuTimes(current)) {
			std::cerr << "\nFailed to get CPU times\n";
			return 1;
		}
		// Вычисляем разницу между измерениями.
		const auto idleDelta =
			current.idle - previous.idle;

		const auto kernelDelta =
			current.kernel - previous.kernel;

		const auto userDelta =
			current.user - previous.user;

		const auto totalDelta =
			kernelDelta + userDelta;

		if (totalDelta > 0 && totalDelta >= idleDelta) {

			const double usage =
				100.0 * (totalDelta - idleDelta) / totalDelta;

			std::cout
				<< "\rCPU Usage: "
				<< usage
				<< "%   "
				<< std::flush;
		}

		// Текущее измерение становится предыдущим.
		previous = current;
	}

	return 0;
}