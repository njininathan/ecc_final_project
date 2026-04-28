#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <time.h>
#include <chrono>
#include <iomanip>

/* Author: Njini Nathan Fofeyin*/

// To be considered - hide constants in header file.
const double rate = 0.5;
const int max_deg = 10;       // Maximum degree of both variable and check nodes
const int max_iter = 20;      // Maximum number of iterations for decoding
const int Z = 96;             // Size of each circulant block
const int H_rows = 1152;      // == num of pairty bits for 1/2 rate code
const int H_cols = 2304;
const int message_len = 1152;
const int parity_len = 1152;
const int blocks = 12;        // Number of 96-bit blocks in parity
const int total_parity = Z * blocks;

// Statistical convergence parameters
const int min_bits_simulated = 10000000;     // Minimum number of information bits to simulate
const int min_errors_required = 100;      // Minimum number of errors before stopping
const int max_num_rounds = 100000;        // Maximum simulation rounds per SNR point
const double confidence_threshold = 0.1;  // Stop when 95% confidence interval is within 10% of mean - tolerance



// Timing functions
std::chrono::high_resolution_clock::time_point tic() {
    return std::chrono::high_resolution_clock::now();
}

double toc(std::chrono::high_resolution_clock::time_point start) {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    return elapsed.count(); // returns total seconds
}

// Utility functions
double randn() {
    static std::random_device rd;      // Used to seed the random number engine
    static std::mt19937 gen(rd());     // Mersenne Twister engine
    static std::normal_distribution<double> dist(0.0, 1.0); // Mean = 0, Stddev = 1

    return dist(gen);
}

float rand49()
{  /*rand_max=7FFF (32767) */
	static int Num=0;
    double number;
    int    i;
    i=rand();
    number=(double)(i)/((unsigned) (RAND_MAX+1));
    Num++;
    if (Num >=RAND_MAX){ // use of conditional unclear. Ask TA.
		time_t t;
		t=time(NULL); // seeder? Verify from TA - it does nothing
        Num=0;
    }
    return (float)number;
}

double Normal()
{
	static int iset=0;
    static double qset;
    double vx,vy,r,temp;
    if (iset==0)//noise=normal*deviate
    {
	    do
        {
		    vx=2.0*rand49()-1.0;
            vy=2.0*rand49()-1.0;
            r =vx*vx+vy*vy;
        }while (r >=1.0 || r==0);
        temp=sqrt(-2.0*log(r)/r);
        qset=vy*temp;
        iset=1;
        return (vx*temp);
    }
    else
    {
   	    iset=0;
        return qset;
    }
}

class Encoder {
public:
    void encode(int* message, int* codeword) {
        std::ifstream fin1("positions_cn_vn.txt");
        if (!fin1) {
            std::cerr << "Failed to open positions_cn_vn.txt\n";
            return;
        }

        int temp[H_rows] = {0}; // store result of H_m*m^T but more efficient using positions
        int parity[parity_len] = {0};

        for (int row = 0; row < H_rows; row++) {
            int bit = 0;
            int col;
            while (fin1 >> col) { // read in the position data - count is from 1 - 2304 correponding to index 0 - 2303 (cols)
                if (col >= H_cols) break; // just in case 
                if (col < message_len) { // this is correct - Please Do Not Touch
                    bit ^= message[col]; // XOR relevant message bit - indexing from 
                }
                // Stop at end of line
                if (fin1.peek() == '\n') break;
            }
            temp[row] = bit; // normally, this is a column vector in mathematical representation
        }

        // Calculate the first set p0 with 96 bits.
        for (int i = 0; i < Z; i++) {
            for (int b = 0; b < blocks; b++) {
                parity[i] ^= temp[b * Z + i]; // block summation equations
            }
        }
        
        // put current message and current parity bit in codeword and utilize this subpart to calculate future parity bits
        for (int i = 0; i < message_len; i++) codeword[i] = message[i];
        for (int i = message_len; i < message_len + Z; i++) codeword[i] = parity[i-message_len];
        
        // double diagonal encoding
        fin1.clear();                 // Clear EOF flag
        fin1.seekg(0, std::ios::beg); // Rewind to beginning

        for (int i = message_len + Z; i < H_cols; i++){ // go through the parity part of the codeword after p0 has been calculated.
            int pbit = 0;
            int col;
            while (fin1 >> col) {
                if (col >= H_cols) break;
                if (col < i) { // up to the highest current known parity bits
                    pbit ^= codeword[col];
                }
                if (fin1.peek() == '\n') break;
            }
            codeword[i] = pbit; // update codeword before you continue to the next parity in the columns of H
        }
        fin1.close();
    }
};

class Channel {
private:
    double sigma; // since only the channel knows
    
public:
    Channel(double sigma_val) : sigma(sigma_val) {}
    
    void addNoise(int* mod_codeword, double* received) { //point to grab these variables
        for (int i = 0; i < H_cols; i++) received[i] = mod_codeword[i] + (sigma*Normal()); 
    }
    
    void modulate(int* codeword, int* mod_codeword) {
        for (int i = 0; i < H_cols; i++) mod_codeword[i] = 1 - 2*codeword[i];
    }
    
    void setSigma(double sigma_val) {
        sigma = sigma_val;
    }
};

class Decoder {
private:
    // Store connectivity information - from file
    int cn_neighbors[H_rows][max_deg];
    int cn_degrees[H_rows];
    int vn_neighbors[H_cols][max_deg];
    int vn_degrees[H_cols];
    double sigma;
    
    void loadConnectivity() {
        std::ifstream fin1("positions_cn_vn.txt");
        std::ifstream fin2("positions_vn_cn.txt");
        
        if (!fin1 || !fin2) {
            std::cerr << "Failed to open connectivity files\n";
            return;
        }
        
        // Initialize degrees
        for (int i = 0; i < H_rows; i++) cn_degrees[i] = 0;
        for (int i = 0; i < H_cols; i++) vn_degrees[i] = 0;
        
        // Initialize neighbor arrays - set -1 to indiciate not connected to store positions from file
        for (int i = 0; i < H_rows; i++) {
            for (int j = 0; j < max_deg; j++) {
                cn_neighbors[i][j] = -1;
            }
        }
        for (int i = 0; i < H_cols; i++) {
            for (int j = 0; j < max_deg; j++) {
                vn_neighbors[i][j] = -1;
            }
        }
        
        // Read CN neighbor information
        for (int i = 0; i < H_rows; i++) {
            int col;
            int degree = 0;
            while (fin1 >> col && degree < max_deg) {
                cn_neighbors[i][degree] = col; // stores locations of 1s for each cn
                degree++;
                if (fin1.peek() == '\n') break;
            }
            cn_degrees[i] = degree;
        }
        fin1.close();

        // Read VN neighbor information
        for (int j = 0; j < H_cols; j++) {
            int row;
            int degree = 0;
            while (fin2 >> row && degree < max_deg) {
                vn_neighbors[j][degree] = row;
                degree++;
                if (fin2.peek() == '\n') break;
            }
            vn_degrees[j] = degree;
        }
        fin2.close();
    }
    
public:
    Decoder(double sigma_val) : sigma(sigma_val) {
        loadConnectivity();
    }
    
    void setSigma(double sigma_val) {
        sigma = sigma_val;
    }
    
    void decode(double* received, int* v) {
        // Initialize
        double L_j[H_cols] = {0};
        for (int i = 0; i < H_cols; i++) L_j[i] = 2*received[i]/pow(sigma,2);
        
        // Message arrays
        double L_ij[H_rows][max_deg] = {0}; // CN to VN messages
        double L_ji[H_cols][max_deg] = {0}; // VN to CN messages

        // Initialize VN to CN messages
        for (int j = 0; j < H_cols; j++) {
            for (int k = 0; k < vn_degrees[j]; k++) {
                L_ji[j][k] = L_j[j];
            }
        }
        
        for (int iter = 0; iter < max_iter; iter++) {
            
        // CN update: Min-Sum version
        for (int i = 0; i < H_rows; i++) {
            for (int k = 0; k < cn_degrees[i]; k++) {
                int vn_idx = cn_neighbors[i][k];
                if (vn_idx == -1) continue;

                double min_abs = 1e9; // initialize to large enough value - this is correct - DNC
                double sign_product = 1.0;

                for (int m = 0; m < cn_degrees[i]; m++) {
                    if (m == k) continue;
                    int other_vn = cn_neighbors[i][m];
                    if (other_vn == -1) continue;

                    // Find L_ji[other_vn][*] going to CN i
                    int other_msg_idx = -1;
                    for (int n = 0; n < vn_degrees[other_vn]; n++) {
                        if (vn_neighbors[other_vn][n] == i) {
                            other_msg_idx = n;
                            break;
                        }
                    }

                    if (other_msg_idx != -1) {
                        double msg = L_ji[other_vn][other_msg_idx];
                        sign_product *= (msg >= 0) ? 1.0 : -1.0;
                        double abs_msg = fabs(msg);
                        if (abs_msg < min_abs) min_abs = abs_msg;
                    }
                }

                // Optional: Offset or normalization to improve performance
                // min_abs = std::max(min_abs - 0.15, 0.0); // offset Min-Sum

                L_ij[i][k] = sign_product * min_abs;
            }
        }

                    
                    // VN update: compute messages from VNs to CNs
                    for (int j = 0; j < H_cols; j++) {
                        // Compute total LLR first
                        double total_llr = L_j[j];
                        for (int k = 0; k < vn_degrees[j]; k++) {
                            int cn_idx = vn_neighbors[j][k];
                            if (cn_idx == -1) break;
                            
                            // Find corresponding message from CN to VN
                            int cn_msg_idx = -1;
                            for (int m = 0; m < cn_degrees[cn_idx]; m++) {
                                if (cn_neighbors[cn_idx][m] == j) {
                                    cn_msg_idx = m;
                                    break;
                                }
                            }
                            
                            if (cn_msg_idx != -1) {
                                total_llr += L_ij[cn_idx][cn_msg_idx];
                            }
                        }
                        
                        // Update VN to CN messages
                        for (int k = 0; k < vn_degrees[j]; k++) {
                            int cn_idx = vn_neighbors[j][k];
                            if (cn_idx == -1) break;
                            
                            // Find corresponding message from CN to VN
                            int cn_msg_idx = -1;
                            for (int m = 0; m < cn_degrees[cn_idx]; m++) {
                                if (cn_neighbors[cn_idx][m] == j) {
                                    cn_msg_idx = m;
                                    break;
                                }
                            }
                            
                            if (cn_msg_idx != -1) {
                                L_ji[j][k] = total_llr - L_ij[cn_idx][cn_msg_idx];
                            } else {
                                L_ji[j][k] = total_llr;
                            }
                        }
                    }
                    
            // LLR total and Make hard decisions
            for (int j = 0; j < H_cols; j++) {
                double total_llr = L_j[j];
                for (int k = 0; k < vn_degrees[j]; k++) {
                    int cn_idx = vn_neighbors[j][k];
                    if (cn_idx == -1) break;
                    
                    // Find corresponding message from CN to VN
                    int cn_msg_idx = -1;
                    for (int m = 0; m < cn_degrees[cn_idx]; m++) {
                        if (cn_neighbors[cn_idx][m] == j) {
                            cn_msg_idx = m;
                            break;
                        }
                    }
                    
                    if (cn_msg_idx != -1) {
                        total_llr += L_ij[cn_idx][cn_msg_idx];
                    }
                }
                
                v[j] = (total_llr < 0) ? 1 : 0;
            }
            
            // Check stopping criterion - reused for checking codeword validity //- remember to convert to function
            bool all_satisfied = true;
            for (int i = 0; i < H_rows; i++) {
                int parity_check = 0;
                for (int k = 0; k < cn_degrees[i]; k++) {
                    int vn_idx = cn_neighbors[i][k];
                    if (vn_idx == -1) break;
                    parity_check ^= v[vn_idx];
                }
                if (parity_check != 0) {
                    all_satisfied = false;
                    break;
                }
            }
            
            if (all_satisfied) {
                break;
            }
        }
    }
};

// Statistics tracking class
class StatisticsTracker { // calculate mean and check for convergence within 95% interval
private:
    std::vector<double> error_rates; //double error_rates[100];
    int window_size;
    
public:
    StatisticsTracker(int window = 100) : window_size(window) {}
    
    void addErrorRate(double rate) {
        error_rates.push_back(rate);
        if (error_rates.size() > window_size) {
            error_rates.erase(error_rates.begin()); // replace the ber value at the beginning with the latest
        }
    }
    
    bool hasConverged() {
        if (error_rates.size() < window_size) return false;
        
        // Calculate mean and standard deviation of recent error rates
        double sum = 0.0;
        for (double rate : error_rates) {
            sum += rate;
        }
        double mean = sum / error_rates.size();
        
        double variance = 0.0;
        for (double rate : error_rates) {
            variance += (rate - mean) * (rate - mean);
        }
        variance /= (error_rates.size() - 1);
        double std_dev = sqrt(variance);
        
        // 95% confidence interval half-width - z value of 1.96 directly from table
        double confidence_interval = 1.96 * std_dev / sqrt(error_rates.size());
        
        // Check if confidence interval is within threshold of mean
        return (confidence_interval / mean) < confidence_threshold;
    }
};

class BERSimulator {
private:
    Encoder encoder;
    Decoder decoder;
    Channel channel;
    double Eb_N0_dB;
    double Eb_N0;
    double sigma;
    
public:
    BERSimulator(double Eb_N0_dB_val) : 
        Eb_N0_dB(Eb_N0_dB_val),
        Eb_N0(pow(10.0, Eb_N0_dB_val/10.0)),
        sigma(sqrt(0.5/(pow(10.0, Eb_N0_dB_val/10.0)*rate))),
        decoder(sqrt(0.5/(pow(10.0, Eb_N0_dB_val/10.0)*rate))),
        channel(sqrt(0.5/(pow(10.0, Eb_N0_dB_val/10.0)*rate))) {}
    
    void setEbN0(double Eb_N0_dB_val) {
        Eb_N0_dB = Eb_N0_dB_val;
        Eb_N0 = pow(10.0, Eb_N0_dB_val/10.0);
        sigma = sqrt(0.5/(Eb_N0*rate));
        decoder.setSigma(sigma);
        channel.setSigma(sigma);
    }
    
    struct SimulationResults { // group result variables
        int total_errors;
        long long bits_simulated; // millions in size
        int rounds_completed;
        double ber;
        bool converged; // convergence checker
        double elapsed_time;
    };
    
    SimulationResults runSimulation() {
        SimulationResults results = {0, 0, 0, 0.0, false, 0.0}; // for the above 6 parameters
        auto start_time = tic();
        
        StatisticsTracker stats_tracker;
        int error_check_interval = 50; // Check convergence every 50 rounds 
        
        std::cout << "Simulating Eb/N0 = " << std::fixed << std::setprecision(1) << Eb_N0_dB << " dB..." << std::endl; // nice display in terminal
        
        for (int round = 0; round < max_num_rounds; round++) {
            // Generate random message
            int message[message_len] = {0};
            for (int i = 0; i < message_len; i++) {
                message[i] = (Normal() < 0.0) ? 1 : 0;
            }

            int codeword[H_cols] = {0};
            
            // Encode
            encoder.encode(message, codeword);
            
            // Modulate
            int mod_codeword[H_cols];
            channel.modulate(codeword, mod_codeword);
            
            // Add channel noise
            double received[H_cols] = {0};
            channel.addNoise(mod_codeword, received);
            
            // Decode
            int decoded[H_cols] = {0};
            decoder.decode(received, decoded);
            
            // Count bit errors (only in information bits for BER calculation)
            int bit_errors = 0;
            for (int i = 0; i < message_len; i++) {
                if (decoded[i] != message[i]) {
                    bit_errors++;
                }
            }
            
            results.total_errors += bit_errors;
            results.bits_simulated += message_len;
            results.rounds_completed = round + 1;
            
            // Check convergence periodically
            if ((round + 1) % error_check_interval == 0) {
                double current_ber = (double)results.total_errors / results.bits_simulated;
                stats_tracker.addErrorRate(current_ber);
                
                // Print progress
                if ((round + 1) % (error_check_interval * 4) == 0) {
                    std::cout << "  Round " << round + 1 << ": " 
                              << results.total_errors << " errors in " 
                              << results.bits_simulated << " bits (BER = " 
                              << std::scientific << std::setprecision(3) << current_ber << ")" << std::endl;
                }
                
                // Check dual stopping criteria
                bool min_bits_reached = results.bits_simulated >= min_bits_simulated;
                bool min_errors_reached = results.total_errors >= min_errors_required;
                bool statistics_converged = stats_tracker.hasConverged();
                
                if (min_bits_reached && min_errors_reached && statistics_converged) {
                    results.converged = true;
                    std::cout << "  --> Statistical convergence achieved!" << std::endl;
                    break;
                }
                
                // For very low BER (high SNR), relax the error requirement if we have enough bits
                if (results.bits_simulated >= min_bits_simulated * 10 && results.total_errors >= 50) { // 100 million points
                    if (statistics_converged) {
                        results.converged = true;
                        std::cout << "  --> Convergence with relaxed error threshold (high SNR)" << std::endl;
                        break;
                    }
                }
            }
        }
        
        results.ber = (double)results.total_errors / results.bits_simulated;
        results.elapsed_time = toc(start_time);
        
        return results;
    }
};

int main() {
    auto total_start = tic();
    srand(42); // Fixed seed for reproducibility
    
    double Eb_N0_dB_values[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2};
    int num_snr_points = sizeof(Eb_N0_dB_values) / sizeof(Eb_N0_dB_values[0]);
    
    std::cout << "LDPC Simulation with Statistical Convergence" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << "Minimum bits per SNR: " << min_bits_simulated << std::endl;
    std::cout << "Minimum errors required: " << min_errors_required << std::endl;
    std::cout << "Maximum rounds per SNR: " << max_num_rounds << std::endl;
    std::cout << std::endl;
    
    // Display results nicely
    // Print heading
    std::cout << std::setw(8) << "Eb/N0(dB)" 
              << std::setw(12) << "Total Err" 
              << std::setw(15) << "Bits Sim"
              << std::setw(10) << "Rounds"
              << std::setw(15) << "BER"
              << std::setw(12) << "Converged"
              << std::setw(10) << "Time(s)" << std::endl;
    std::cout << std::string(82, '-') << std::endl;
    
    for (int i = 0; i < num_snr_points; i++) {
        BERSimulator simulator(Eb_N0_dB_values[i]);
        BERSimulator::SimulationResults results = simulator.runSimulation();
        
        // Print results
        std::cout << std::setw(8) << std::fixed << std::setprecision(1) << Eb_N0_dB_values[i]
                  << std::setw(12) << results.total_errors
                  << std::setw(15) << results.bits_simulated
                  << std::setw(10) << results.rounds_completed
                  << std::setw(15) << std::scientific << std::setprecision(3) << results.ber
                  << std::setw(12) << (results.converged ? "Yes" : "No")
                  << std::setw(10) << std::fixed << std::setprecision(1) << results.elapsed_time << std::endl;
    }
    
    double total_elapsed = toc(total_start);
    std::cout << std::string(82, '=') << std::endl;
    std::cout << "Total simulation time: " << std::fixed << std::setprecision(1) 
              << total_elapsed << " seconds" << std::endl;
    
    return 0;
}