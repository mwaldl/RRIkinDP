#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <numeric> // For accumulate
#include <Eigen/Dense>
#include <boost/program_options.hpp>
#include <sstream>
#include <iterator>

namespace po = boost::program_options;

// Function to read a vector from a file (for the initial state)
Eigen::VectorXd readVector(const std::string &filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<double> values;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::istringstream stream(line);
            double value;
            stream >> value;
            values.push_back(value);
        }
    } else {
        throw std::runtime_error("Could not open file: " + filename);
    }

    Eigen::VectorXd vector(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        vector(i) = values[i];
    }

    return vector;
}

// Function to read rate matrix from a file
Eigen::MatrixXd readMatrix(const std::string &filename) {
    std::ifstream file(filename);
    std::string line;
    std::vector<std::vector<double>> values;
    int rows = 0;
    int cols = 0;

    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::istringstream stream(line);
            std::vector<double> row_values;
            double value;

            while (stream >> value) {
                row_values.push_back(value);
            }

            if (rows == 0) {
                cols = row_values.size();
            }
            values.push_back(row_values);
            rows++;
        }
    } else {
        throw std::runtime_error("Could not open file: " + filename);
    }

    Eigen::MatrixXd matrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix(i, j) = values[i][j];
        }
    }
    // matrix.transposeInPlace();  // Adjust for transpose if needed

    return matrix;
}

class MarkovChainSimulator {
public:
    MarkovChainSimulator(const Eigen::MatrixXd& rateMatrix, Eigen::VectorXd& initialDistribution,
                         double startTime, double endTime, double factor, const std::string& outputFile)
        : R(rateMatrix),
          initialDistribution(initialDistribution),
          startTime(startTime), endTime(endTime), factor(factor), outputFile(outputFile) {
        symmetrizeMatrix();
    }

    void simulate() {
        std::vector<double> times;
        double currentTime = startTime;

        // Generate time points using the exponential factor
        while (currentTime <= endTime) {
            times.push_back(currentTime);
            currentTime *= factor;  // Update to the next time point
        }

        std::ofstream outFile(outputFile);
        outFile << std::fixed << std::setprecision(6);
        outFile << "Time";
        for (int i = 0; i < initialDistribution.size(); ++i) {
            outFile << " P(" << i + 1 << ")";
        }
        outFile << "\n";

        for (double t : times) {
            Eigen::VectorXd probabilities = computeProbabilities(t);
            outFile << t;
            for (int i = 0; i < probabilities.size(); ++i) {
                outFile << " " << probabilities[i];
            }
            outFile << "\n";
        }

        outFile.close();
    }

    // Function to get the equilibrium distribution
    Eigen::VectorXd getEquilibriumDistribution() const {
        return computeEquilibriumDistribution();
    }

private:
    Eigen::MatrixXd R;  // Rate matrix
    Eigen::MatrixXd U;  // Symmetrized matrix
    Eigen::VectorXd initialDistribution;
    double startTime;
    double endTime;
    double factor;  // Factor for exponential time stepping
    std::string outputFile;

    void symmetrizeMatrix() {

    }

    Eigen::VectorXd computeEquilibriumDistribution() const {

    }

    Eigen::VectorXd computeProbabilities(double time) {

    }
};

int main(int argc, char* argv[]) {
    std::string ratesFile;
    std::string statesFile;
    std::string outputFile;
    double startTime = 0.1;
    double endTime = 10.0;
    double factor = 1.02;  // Default factor for exponential growth

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Help screen")
        ("ratesfile,r", po::value<std::string>(&ratesFile)->required(), "Path to the rate matrix file")
        ("states,s", po::value<std::string>(&statesFile)->required(), "Path to the initial state distribution file")
        ("output,o", po::value<std::string>(&outputFile)->default_value("output.txt"), "Output file")
        ("starttime,t", po::value<double>(&startTime)->default_value(0.1), "Start time for simulation")
        ("endtime,e", po::value<double>(&endTime)->default_value(10.0), "End time for simulation")
        ("factor,f", po::value<double>(&factor)->default_value(1.02), "Exponential growth factor for time stepping");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    try {
        Eigen::VectorXd initialDistribution = readVector(statesFile);
        Eigen::MatrixXd rateMatrix = readMatrix(ratesFile);

        MarkovChainSimulator simulator(rateMatrix, initialDistribution, startTime, endTime, factor, outputFile);
        simulator.simulate();

        Eigen::VectorXd equilibriumDistribution = simulator.getEquilibriumDistribution();
        std::cout << "Equilibrium Distribution: " << equilibriumDistribution.transpose() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}


