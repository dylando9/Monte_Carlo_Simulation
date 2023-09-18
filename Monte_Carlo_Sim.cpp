#include <iostream>
#include <fstream> // for reading/writing files for read_xyz function
#include <array>   // for std::array for read_xyz function
#include <vector>  // for std::vector for read_xyz function
#include <utility> // for std::pair for read_xyz function
#include <cmath>   // for the calculate_tail_correction
#include <random>  // for random numbers for the random_double function
#include <chrono>  // for generating the seed for the random_double function

std::default_random_engine re; // Random Double Function

/*! Generate a random double within a given range */
double random_double(const double &lower_bound, const double &upper_bound)
{
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    return dist(re);
}

/*! Generate a random integer within a given range
    The generated integer will be on the range [a,b)
*/
double random_integer(const int &lower_bound, const int &upper_bound)
{
    // dist will return [a,b] but we want [a,b)
    std::uniform_int_distribution<int> dist(lower_bound, upper_bound - 1);
    return dist(re);
}

typedef std::array<double, 3> AtomCoord;    // Make array with 3 columns
typedef std::vector<AtomCoord> Coordinates; // Make the vector of the array

// Configure random Points
double configure_random_points(const int &number_of_atoms, const double &density)
{
    double box_volume = number_of_atoms / density;
    double box_length = pow(box_volume, (1 / 3));

    // std::vector<double> total_coordinates;
    Coordinates coords;

    // Generate the x, y, z coordinates based on the number of atoms given and put them in the total coordinates vector
    for (int i = 0; i < number_of_atoms; i++)
    {
        double x = random_double(0.0, 1.0);
        double y = random_double(0.0, 1.0);
        double z = random_double(0.0, 1.0);
        AtomCoord coord;

        // coord[0] >> coord[1] >> coord[2];
        coord[0] = x;
        coord[1] = y;
        coord[2] = z;

        coords.push_back(coord);
    }
    return coords, box_length;
}

////////////////////////////////////////// calculate_distance FUNCTION /////////////////////////////////////////////////////////////////////

double calculate_distance(const AtomCoord &coord1, const AtomCoord &coord2, double box_length = 0.0)
{

    /*
    This first part will replace the dim_dist = coord1[i] - coord2[i] we had in python
    if the box length is not provided as input for the function then we will just calculate
    the distance like normal.
    */

    double dx = coord1[0] - coord2[0];
    double dy = coord1[1] - coord2[1];
    double dz = coord1[2] - coord2[2];

    /* However if the box length is provided and is not 0
    we will take into account the periodic boundaries using this if statement.
    */

    if (box_length > 0.0)
    {
        dx = dx - box_length * std::round(dx / box_length);
        dy = dy - box_length * std::round(dy / box_length);
        dz = dz - box_length * std::round(dz / box_length);
    }

    double distance = std::sqrt(dx * dx + dy * dy + dz * dz); // combines the ^2 in sqrt in one line vs 2 in ython
    return distance;
}

///////////////////////////////////////////////// calculate_LJ FUNCTION ////////////////////////////

double calculate_LJ(const double &r_ij)
{

    /*
    The LJ interaction energy between two particles.

        Computes the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

        Parameters
        ----------
        r_ij : double
            The distance between the particles in reduced units.

        Returns
        -------
        pairwise_energy : double
            The pairwise Lennard Jones interaction energy in reduced units.

    */

    double r6_term = (1 / r_ij) * (1 / r_ij) * (1 / r_ij) * (1 / r_ij) * (1 / r_ij) * (1 / r_ij);
    double r12_term = (r6_term) * (r6_term);
    double pairwise_energy = 4.0 * (r12_term - r6_term);

    return pairwise_energy;
}

///////////////////////////////////////////////// read_xyz FUNCTION ////////////////////////////////////////////////////////////////////

std::pair<Coordinates, double> read_xyz(const std::string &file_path)
{
    // Opens up a file stream for input
    std::ifstream infile(file_path);

    // Check that it was successfully opened
    if (!infile.is_open())
    {
        throw std::runtime_error("File path in read_xyz does not exist!");
    }

    double dummy; // Data that is thrown away (box length, atom indices)
    double box_length;
    int num_atoms;

    // Grab box_length from first number, throw the rest away
    infile >> box_length >> dummy >> dummy;

    // now the number of atoms
    infile >> num_atoms;

    // Uncomment to help troubleshoot
    // std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;

    // Stores the atomic coordinates
    // Remember, this is a vector of arrays
    Coordinates coords;

    for (int i = 0; i < num_atoms; i++)
    {
        AtomCoord coord;

        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];

        // Add to the vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}

/////////////////////////////////////////////accept_reject FUNCTION//////////////////////////////////////////////////////////////////////////////

bool accept_or_reject(const double &delta_e, const double &beta)
{
    if (delta_e <= 0)
    {
        return true;
    }
    else
    {
        re.seed(std::chrono::system_clock::now().time_since_epoch().count());

        double test = random_double(0.0, 1.0); // THIS IS HOW YOU CALL THE FUNCTION
        // std::cout << test << std::endl;
        double p_acc = std::exp(-beta * delta_e);
        // std::cout << p_acc << std::endl;
        if (test < p_acc)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

/////////////////////////////////////////////calculate_tail_correction FUNCTION///////////////////////////////////////////////

double calculate_tail_correction(const int &num_particles, const double &box_length, const double &cutoff)
{
    double const1 = (8 * M_PI * pow(num_particles, 2)) / (3 * pow(box_length, 3));
    double const2 = ((1.0 / 3.0) * pow((1 / cutoff), 9)) - pow((1 / cutoff), 3);
    return const1 * const2;
}

///////////////////////////////////////////calculate_total_energy FUNCTION///////////////////////////////////////////////////

double calculate_total_energy(const Coordinates &coordinates, const double &box_length, const double &cutoff)
{

    int num_atoms = coordinates.size();
    double total_energy = 0.0; // defining this here so it works in the loop

    for (int i = 0; i < num_atoms; i++)
    {
        for (int j = i + 1; j < num_atoms; j++)
        {
            double dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length);
            {
                if (dist_ij < cutoff)
                {
                    double interaction_energy = calculate_LJ(dist_ij);
                    total_energy += interaction_energy;
                }
            }
        }
    }

    return total_energy;
}

///////////////////////////////////////////////////calculate_pair_energy FUNCTION ////////////////////////////////////////////

double calculate_pair_energy(const Coordinates &Coords, const int &i_particle, const double &box_length, const double &cutoff)
{
    double e_total = 0.0;
    AtomCoord i_position = Coords[i_particle];
    size_t num_atoms = Coords.size();

    for (int j_particle = 0; j_particle < num_atoms; j_particle++)
    {
        if (i_particle != j_particle)
        {
            AtomCoord j_position = Coords[j_particle];
            double rij = calculate_distance(i_position, j_position, box_length);
            if (rij < cutoff)
            {
                double e_pair = calculate_LJ(rij);
                e_total += e_pair;
            }
        }
    }
    return e_total;
}

/////////////////////////////////////////////////run_simulation FUNCTION//////////////////////////////////////////////////////////////////////

std::vector<Coordinates> run_simulation(Coordinates coordinates, const double &box_length, const double &cutoff, const double &reduced_temperature, const int &num_steps, const double &max_displacement, const int &freq = 1000)
{
    std::vector<int> steps;
    std::vector<double> energies;
    std::vector<Coordinates> all_coordinates;

    double beta = 1.0 / reduced_temperature;

    size_t num_particles = coordinates.size();

    double total_energy = calculate_total_energy(coordinates, box_length, cutoff);
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff);

    std::cout << "The starting energy is " << total_energy << " ." << std::endl;

    for (int step = 0; step < num_steps; step++)
    {
        int random_particle = random_double(0, num_particles);

        double current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);

        coordinates[random_particle][0] += x_rand;
        coordinates[random_particle][1] += y_rand;
        coordinates[random_particle][2] += z_rand;

        double proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double delta_energy = (proposed_energy - current_energy);

        bool accept = accept_or_reject(delta_energy, beta);

        if (accept == true)
        {
            total_energy += delta_energy;
        }
        else
        {
            coordinates[random_particle][0] -= x_rand;
            coordinates[random_particle][1] -= y_rand;
            coordinates[random_particle][2] -= z_rand;
        }

        if (step % freq == 0)
        {
            std::cout << step << "  " << total_energy / num_particles << std::endl;
            steps.push_back(step);
            energies.push_back(total_energy / num_particles);
            all_coordinates.push_back(coordinates);
        }
    }

    return all_coordinates;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(void)
{

    double reduced_temperature = 1.5;
    int num_steps = 50000;
    double max_displacement = 0.1;
    double cutoff = 3.0;
    double freq = 1000.0;

    std::string file_path = "sample_config1.txt";

    auto p = read_xyz(file_path);

    std::cout << p.first[0][0] << " : " << p.second << std::endl;

    auto all_coordinates = run_simulation(p.first, p.second, cutoff, reduced_temperature, num_steps, max_displacement, freq = 1000);

    return 0;
}