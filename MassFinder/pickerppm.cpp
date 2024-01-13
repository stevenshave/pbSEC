#include <cctype>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "Usage: ./picker resultsfile wellfilename tolerance(in ppm) mass1 "
                 "mass2 massn....\n";
    exit(-1);
  }

  std::string resultsfile = argv[1];
  std::string wellFilename = argv[2];

  float tol = atof(argv[3]);


  std::vector<float> masses, counts;

  for (int i = 4; i < argc; i++) {
    masses.push_back(atof(argv[i]));
    counts.push_back(0.0);
  }
  std::vector<std::vector<float> > massFrequency (100, std::vector<float>(masses.size(), 0.0f));

  std::cerr << "Tol=" << tol << " ppm\n";

  std::string line;
  size_t pos1, pos2;
  unsigned int scannumber;
  std::string mslevel;
  float mass;
  float intensity;

  unsigned maxtime = 120;


  std::getline(std::cin, line);
  while (std::cin && !std::cin.eof()) {
    std::getline(std::cin, line);
    if (line.find("ms1 ") == std::string::npos)
      continue;

    pos1 = line.find_first_not_of(' ', 0);
    pos2 = line.find_first_of(' ', pos1 + 1);
    scannumber = std::stoi(line.substr(pos1, pos2 - pos1), 0, 10);
    if (scannumber < 120)
      continue;

    pos1 = line.find_first_not_of(' ', pos2);
    pos2 = line.find_first_of(' ', pos1 + 1);

    pos1 = line.find_first_not_of(' ', pos2);
    pos2 = line.find_first_of(' ', pos1 + 1);
    mass = std::stof(line.substr(pos1, pos2 - pos1), 0);

    pos1 = line.find_first_not_of(' ', pos2);
    pos2 = line.find_first_of(' ', pos1 + 1);
    intensity = std::stof(line.substr(pos1, pos2 - pos1), 0);
    if (intensity == 0)
      continue;

	float tolForMass;
    for (size_t i = 0; i < masses.size(); ++i) {
		tolForMass = (masses[i] * tol) / static_cast<float>(1e6);
      if ((mass >= (masses[i] - tolForMass)) && (mass <= (masses[i] + tolForMass))) {
        counts[i] += intensity;
		maxtime = std::max(maxtime, scannumber);
		massFrequency[(scannumber - 120) / 60][i] += intensity;
      }
    }
  }

  std::ofstream outfile;
  outfile.open(resultsfile, std::ios_base::app);
  for (int i = 0; i < masses.size(); ++i) {
	  outfile << wellFilename << "," << masses[i] << "," << counts[i]<<",";
	  for (int j = 0; j < massFrequency.size() && j<=(maxtime-120)/60; ++j) {
		  outfile << massFrequency[j][i] << ",";
	  }
	  outfile << "\n";
  }

  return 0;
}
