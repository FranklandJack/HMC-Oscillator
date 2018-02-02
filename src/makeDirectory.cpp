#include "makeDirectory.hpp"

std::string makeDirectory()
{
	// Begin a timer for the output directory.
    auto start = std::chrono::system_clock::now();

	// Create string from the time the program started.
	time_t startTime = std::chrono::system_clock::to_time_t(start);
	std::string outputName = ctime(&startTime);

	// Strip out and replace difficult characters.
	std::transform(outputName.begin(), outputName.end(), outputName.begin(), [](char ch) {return ch == ' ' ? '_' : ch;});
	std::transform(outputName.begin(), outputName.end(), outputName.begin(), [](char ch) {return ch == ':' ? '-' : ch;});
	outputName.erase(std::remove(outputName.begin(), outputName.end(), '\n'), outputName.end());

	// Create directory path from the string.
	boost::filesystem::path outPath = outputName;

	// If user calls program more than once a second so that directories will be overwritten append an index.
	for(int i = 2; boost::filesystem::exists(outPath) && i < 100; ++i)
	{
	    std::stringstream ss;
	    ss << outPath << "(" << i << ")";
	    outPath = ss.str();
	}

	// Create the directory for output.
	boost::filesystem::create_directory(outPath);

	// Return the string so that the user can create files in the directory.
	return outputName;

}