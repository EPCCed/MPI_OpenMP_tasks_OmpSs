#include <cassert>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sys/time.h>

#include "common/matrix.hpp"
#include "common/heat.hpp"

int initialize(HeatConfiguration &conf, int rowBlocks, int colBlocks, int rowBlockOffset)
{
	conf.matrix = (block_t *) malloc(rowBlocks * colBlocks * sizeof(block_t));
	if (conf.matrix == NULL) {
		fprintf(stderr, "Memory cannot be allocated!\n");
		exit(1);
	}
	
	initializeMatrix(conf, conf.matrix, rowBlocks, colBlocks, rowBlockOffset);
	
	return 0;
}

int finalize(HeatConfiguration &conf)
{
	assert(conf.matrix != nullptr);
	free(conf.matrix);
	conf.matrix = nullptr;
	
	return 0;
}

int writeImage(std::string imageFileName, block_t *matrix, int rowBlocks, int colBlocks)
{
	// RGB table
	unsigned int r[1024], g[1024], b[1024];
	
	// Prepare the RGB table
	int n = 1023;
	for (int i = 0; i < 256; i++) {
		r[n] = 255; g[n] = i; b[n] = 0;
		n--;
	}
	
	for (int i = 0; i < 256; i++) {
		r[n] = 255 - i; g[n] = 255; b[n] = 0;
		n--;
	}
	
	for (int i = 0; i < 256; i++) {
		r[n] = 0; g[n] = 255; b[n] = i;
		n--;
	}
	
	for (int i = 0; i < 256; i++) {
		r[n] = 0; g[n] = 255 - i; b[n] = 255;
		n--;
	}
	
	// Find minimum and maximum
	double min = DBL_MAX;
	double max = -DBL_MAX;
	traverseByRows(matrix, rowBlocks, colBlocks,
		[&](int x, int y, double value) {
			if (value > max)
				max = value;
			if (value < min)
				min = value;
		}
	);
	
	int rows = (rowBlocks - 2) * BSX + 2;
	int cols = (colBlocks - 2) * BSY + 2;
	std::ofstream file;
	file.open(imageFileName);
	
	file << "P3" << std::endl;
	file << cols << " " << rows << std::endl;
	file << 255 << std::endl;
	
	traverseByRows(matrix, rowBlocks, colBlocks,
		[&](int x, int y, double value) {
			int k = 0;
			if (max - min != 0) {
				k = (int)(1023.0 * (value - min)/(max - min));
			}
			file << r[k] << " " << g[k] << " " << b[k] << "  ";
			if (y == cols - 1) file << std::endl;
		}
	);
	
	file.close();
	
	return 0;
}

void readParameters(int argc, char **argv, HeatConfiguration &conf)
{
	static struct option long_options[] = {
		{"size",         required_argument,  0, 's'},
		{"rows",         required_argument,  0, 'r'},
		{"cols",         required_argument,  0, 'c'},
		{"timesteps",    required_argument,  0, 't'},
		{"sources-file", required_argument,  0, 'f'},
		{"output",       optional_argument,  0, 'o'},
		{"help",         no_argument,        0, 'h'},
		{0, 0, 0, 0}
	};
	
	int c;
	int index;
	while ((c = getopt_long(argc, argv, "ho::fs:r:c:t:", long_options, &index)) != -1) {
		switch (c) {
			case 'h':
				fprintf(stdout, "Usage: %s [OPTION]...\n", argv[0]);
				fprintf(stdout, "Available options:\n");
				fprintf(stdout, "  -s, --size=SIZE\t\tuse SIZExSIZE matrix as the surface (default: 2048)\n");
				fprintf(stdout, "  -r, --rows=ROWS\t\tuse ROWS as the number of rows of the surface (default: 2048)\n");
				fprintf(stdout, "  -c, --cols=COLS\t\tuse COLS as the number of columns of the surface (default: 2048)\n");
				fprintf(stdout, "  -t, --timesteps=TIMESTEPS\tuse TIMESTEPS as the number of timesteps (default: 100)\n");
				fprintf(stdout, "  -f, --sources-file=NAME\tget the heat sources from the NAME configuration file (default: heat.conf)\n");
				fprintf(stdout, "  -o, --output[=NAME]\t\tsave the computed matrix to a PPM file, being 'heat.ppm' the default name (disabled by default)\n");
				fprintf(stdout, "  -h, --help\t\t\tdisplay this help and exit\n");
				exit(0);
			case 'f':
				conf.confFileName = optarg;
				break;
			case 'o':
				conf.generateImage = true;
				if (optarg) {
					conf.imageFileName = optarg;
				}
				break;
			case 's':
				conf.rows = atoi(optarg);
				conf.cols = atoi(optarg);
				break;
			case 'r':
				conf.rows = atoi(optarg);
				break;
			case 'c':
				conf.cols = atoi(optarg);
				break;
			case 't':
				conf.timesteps = atoi(optarg);
				break;
			case '?':
				exit(1);
			default:
				abort();
		}
	}
}

HeatConfiguration readConfiguration(int argc, char **argv)
{
	// Default configuration
	HeatConfiguration conf;
	
	// Read the execution parameters
	readParameters(argc, argv, conf);
	
	std::string line;
	std::ifstream file(conf.confFileName);
	if (!file.is_open()) {
		fprintf(stderr, "Configuration file %s not found!\n", conf.confFileName.c_str());
		exit(1);
	}
	
	std::getline(file, line);
	int n = std::sscanf(line.c_str(), "%d", &(conf.numHeatSources));
	if (n != 1) {
		fprintf(stderr, "Configuration file not correct!\n");
		exit(1);
	}
	
	conf.heatSources = (HeatSource *) malloc(sizeof(HeatSource) * conf.numHeatSources);
	assert(conf.heatSources != nullptr);
	
	for (int i = 0; i < conf.numHeatSources; i++) {
		std::getline(file, line);
		n = std::sscanf(line.c_str(), "%f %f %f %f",
			&(conf.heatSources[i].row),
			&(conf.heatSources[i].col),
			&(conf.heatSources[i].range),
			&(conf.heatSources[i].temperature));
		
		if (n != 4) {
			fprintf(stderr, "Configuration file not correct!\n");
			exit(1);
		}
	}
	
	file.close();
	
	return conf;
}

void refineConfiguration(HeatConfiguration &conf, int rowValue, int colValue)
{
	assert(conf.rows > 0);
	assert(conf.cols > 0);
	assert(conf.timesteps > 0);
	
	if (conf.rows % rowValue) {
		// Make the number of rows divisible by the value
		fprintf(stderr, "Warning: The number of rows (%d) is not divisible by %d. Rounding it...\n", conf.rows, rowValue);
		conf.rows = round(conf.rows, rowValue);
	}
	if (conf.cols % colValue) {
		// Make the number of cols divisible by the value
		fprintf(stderr, "Warning: The number of cols (%d) is not divisible by %d. Rounding it...\n", conf.cols, colValue);
		conf.cols = round(conf.cols, colValue);
	}
}

void printConfiguration(const HeatConfiguration &conf)
{
	fprintf(stdout, "Rows x Cols       : %u x %u\n", conf.rows, conf.cols);
	fprintf(stdout, "Timesteps         : %u\n", conf.timesteps);
	fprintf(stdout, "Num. heat sources : %u\n", conf.numHeatSources);
	
	for (int i = 0; i < conf.numHeatSources; i++) {
		fprintf(stdout, "  %2d: (%2.2f, %2.2f) %2.2f %2.2f \n", i+1,
			conf.heatSources[i].row,
			conf.heatSources[i].col,
			conf.heatSources[i].range,
			conf.heatSources[i].temperature
		);
	}
}

void initializeMatrix(const HeatConfiguration &conf, block_t *matrix, int rowBlocks, int colBlocks, int rowBlockOffset)
{
	const int totalRowBlocks = conf.rowBlocks + 2;
	const int totalRows = conf.rows + 2;
	const int numRows = (rowBlocks - 2) * BSX + 2;
	const int numCols = (colBlocks - 2) * BSY + 2;
	const int rowOffset = rowBlockOffset * BSX;
	
	// Set all elements to zero
	traverseByRows(matrix, rowBlocks, colBlocks,
		[&](int x, int y, double &value) {
			value = 0;
		}
	);
	
	for (int i = 0; i < conf.numHeatSources; i++) {
		const HeatSource &src = conf.heatSources[i];
		
		// Initialize top row
		if (rowBlockOffset == 0) {
			traverseRow(matrix, rowBlocks, colBlocks, 0, 0, numCols,
				[&](int x, int y, double &value) {
					double dist = sqrt(pow((double)y / (double)numCols - src.col, 2) + pow(src.row, 2));
					if (dist <= src.range) {
						value += (src.range - dist) / src.range * src.temperature;
					}
				}
			);
		}
		
		// Initialize bottom row
		if (rowBlockOffset + rowBlocks == totalRowBlocks) {
			traverseRow(matrix, rowBlocks, colBlocks, numRows - 1, 0, numCols,
				[&](int x, int y, double &value) {
					double dist = sqrt(pow((double)y / (double)numCols - src.col, 2) + pow(1 - src.row, 2));
					if (dist <= src.range) {
						value += (src.range - dist) / src.range * src.temperature;
					}
				}
			);
		}
		
		// Initialize left column
		traverseCol(matrix, rowBlocks, colBlocks, 0, 1, numRows - 1,
			[&](int x, int y, double &value) {
				double dist = sqrt(pow(src.col, 2) + pow((double)(rowOffset + x)/(double)totalRows - src.row, 2));
				if (dist <= src.range) {
					value += (src.range - dist) / src.range * src.temperature;
				}
			}
		);
		
		// Initialize right column
		traverseCol(matrix, rowBlocks, colBlocks, numCols - 1, 1, numRows - 1,
			[&](int x, int y, double &value) {
				double dist = sqrt(pow(1 - src.col, 2) + pow((double)(rowOffset + x)/(double)totalRows - src.row, 2));
				if (dist <= src.range) {
					value += (src.range - dist) / src.range * src.temperature;
				}
			}
		);
	}
}

double get_time()
{
	struct timeval tv;
	gettimeofday(&tv, 0);
	
	return tv.tv_sec + 1e-6 * tv.tv_usec;
}

