#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void read_positions(double *positions, int N) {
    FILE *file = fopen("particles_1024.txt", "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening particles.txt\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; ++i) {
        int index;
        double x, y, z;
        if (fscanf(file, "%d %lf %lf %lf", &index, &x, &y, &z) != 4) {
            fprintf(stderr, "Error reading particle coordinates from particles.txt\n");
            exit(EXIT_FAILURE);
        }

        // Assuming positions array is laid out as [x0, y0, z0, x1, y1, z1, ..., xn, yn, zn]
        positions[i * 3] = x;
        positions[i * 3 + 1] = y;
        positions[i * 3 + 2] = z;
        }
        
        fclose(file);
}

void print_forces(double *forces, int N) {
    printf("Forces array:\n");

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("Forces between particles %d and %d:\n", i, j);

            for (int k = 0; k < 3; ++k) {
                printf("   Force[%d][%d] in direction %c: %lf\n", i, j, 'x' + k, forces[i * N * 3 + j * 3 + k]);
            }
        }
    }
}

void print_forces_single(double *forces_single, int N) {
    printf("Forces_single array:\n");

    for (int i = 0; i < N; ++i) {
        printf("Forces on particle %d:\n", i);

        for (int direction = 0; direction < 3; ++direction) {
            printf("   Force[%d] in direction %c: %lf\n", i, 'x' + direction, forces_single[i * 3 + direction]);
        }
    }
}

void calculate_distances(double *positions, double *distances, int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dx = positions[i * 3] - positions[j * 3];
            double dy = positions[i * 3 + 1] - positions[j * 3 + 1];
            double dz = positions[i * 3 + 2] - positions[j * 3 + 2];

            double distance = sqrt(dx * dx + dy * dy + dz * dz);

            // Fill the symmetric matrix
            distances[i * N + j] = distance;
            distances[j * N + i] = distance;
        }
      
        // Diagonal elements are 0
        distances[i * N + i] = 0.0;
    }
}

void calculate_potentials(double *distances, double *potentials, int N) {
    const double r_c = 10.0;

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double distance = distances[i * N + j];

            // Check if distance is within the cutoff radius
            if (distance <= r_c) {
                double r_6 = pow(distance, 6);
                double r_12 = r_6 * r_6;
                potentials[i * N + j] = 1.0 / r_12 - 2.0 / r_6;
                potentials[j * N + i] = potentials[i * N + j]; // Due to symmetry
            } else {
                potentials[i * N + j] = 0.0;
                potentials[j * N + i] = 0.0; // Due to symmetry
            }
        }
       
        // Diagonal elements are 0
        potentials[i * N + i] = 0.0;
    }
}

double calculate_potentials_single(double *distances, int N, int index_i, int index_j) {
    const double r_c = 10.0;
    double potential = 0.0;

    double distance = distances[index_i * N + index_j];

    // Check if distance is within the cutoff radius
    if (distance <= r_c) {
        double r_6 = pow(distance, 6);
        double r_12 = r_6 * r_6;
        potential = 1.0 / r_12 - 2.0 / r_6;
    } else {
        potential = 0.0;
    }
       
    // Diagonal elements are 0
    if(index_i == index_j){
        potential = 0.0;
    }
    return potential;
}

//fills the total_energy array for all particles
//requires the potentials between all particles as input
void calculate_tot_energy(double *potentials, int N, double *total_energy) {
    for (int i = 0; i < N; ++i) {
        total_energy[i] = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                total_energy[i] += potentials[i * N + j];
            }
        }
        total_energy[i] *= 0.5;  // Multiply by 0.5 as per the formula
    }
}                                                                                                                                   

//this function calculates potentials on the fly
//fills the total_energy array for all particles
void calc_tot_energy(double *positions, double *distances, int N, double *total_energy) {
    const double r_c = 10.0;

    // Allocate memory for potentials
    double *potentials = (double *)malloc(N * sizeof(double));
    if (potentials == NULL) {
        fprintf(stderr, "Error allocating memory\n");
        exit(EXIT_FAILURE);
    }
    
    // Calculate potentials and total energy
    for (int i = 0; i < N; ++i) {
        potentials[i] = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                double distance = distances[i * N + j];
    
                // Check if distance is within the cutoff radius
                if (distance <= r_c) {
                    double r_12 = pow(distance, -12);
                    double r_6 = pow(distance, -6);
                    potentials[i] += r_12 - 2.0 * r_6;
                }
            }
        }
    
        total_energy[i] = 0.5 * potentials[i];
    }
    
    // Free allocated memory
    free(potentials);
}

//calculates the total energy for one particle
//and calculates the potential on the fly
double calc_tot_energy_single(double *positions, double *distances, int N, int particle_index) {
    const double r_c = 10.0;

    double potential = 0.0;

    for (int j = 0; j < N; ++j) {
        if (particle_index != j) {
            double distance = distances[particle_index * N + j];

            // Check if distance is within the cutoff radius
            if (distance <= r_c) {
                double r_12 = pow(distance, -12);
                double r_6 = pow(distance, -6);
                potential += r_12 - 2.0 * r_6;
            }
        }
    }
            
    return 0.5 * potential;
}

//this function is wrong
void calc_forces_single(double *positions, double *distances, double *tot_energies, int N, int particle_index, double *forces, double h) {
    const double epsilon = 1e-12; // Default value for h

    // Calculate total energy for the current coordinates
    double energy_current = tot_energies[particle_index];
    
    // Perturb the position forward in x and calculate the total energy
    positions[particle_index * 3] += h;
    double energy_forward_x = calc_tot_energy_single(positions, distances, N, particle_index);
    forces[particle_index * 3] = -(energy_forward_x - energy_current) / h;
    positions[particle_index * 3] -= h; // Reset x position
    
    // Perturb the position forward in y and calculate the total energy
    positions[particle_index * 3 + 1] += h;
    double energy_forward_y = calc_tot_energy_single(positions, distances, N, particle_index);
    forces[particle_index * 3 + 1] = -(energy_forward_y - energy_current) / h;
    positions[particle_index * 3 + 1] -= h; // Reset y position
    
    // Perturb the position forward in z and calculate the total energy
    positions[particle_index * 3 + 2] += h;
    double energy_forward_z = calc_tot_energy_single(positions, distances, N, particle_index);
    forces[particle_index * 3 + 2] = -(energy_forward_z - energy_current) / h;
    positions[particle_index * 3 + 2] -= h; // Reset z position
}

double calc_force_single_direction(double *positions, double *distances, int N, int index_i, int index_j, int direction) {
    const double epsilon = 1e-12;
    double cutoff = 10;
    if(index_i == index_j){
        return 0.0;
    }
    else{
	double distance = distances[index_i * N + index_j];
	if(distance < cutoff){
            // Calculate potential energy for the original configuration
            double potential_original = calculate_potentials_single(distances, N, index_i, index_j);
        
            // Perturb the position of particle i in the specified direction
            positions[index_i * 3 + direction] += epsilon;
        
            // Calculate the new distance to particle j
            double dx = positions[index_i * 3] - positions[index_j * 3];
            double dy = positions[index_i * 3 + 1] - positions[index_j * 3 + 1];
            double dz = positions[index_i * 3 + 2] - positions[index_j * 3 + 2];
            double new_distance = sqrt(dx * dx + dy * dy + dz * dz);
        
            // Calculate the new potential energy
            double potential_new = 1.0 / pow(new_distance, 12) - 2.0 / pow(new_distance, 6);
        
            // Calculate the force using the forward difference formula
            double force = (potential_new - potential_original) / epsilon;
        
            // Reset the perturbed position
            positions[index_i * 3 + direction] -= epsilon;
        
            return (-1 * force);
	}
	else{
            return 0.0;
	}
    }
}

double calc_force_per_particle(double *positions, double *distances, int N, int index_i, int direction) {
    double total_force = 0.0;

    for (int j = 0; j < N; ++j) {
        // Calculate force for the specified direction using calc_force_single_direction
        double force = calc_force_single_direction(positions, distances, N, index_i, j, direction);

        // Accumulate the forces for each particle
        total_force += force;
    }

    return total_force;
}   

int main(int argc, char *argv[])
{
  int rank, size;
  int N; //Number of particles
  double start_time, end_time;
  double start_time2, end_time2;
  double r_c;
  double h; //step size for forward difference approximation
  double *positions;
  double *distances;
  double *potentials;
  double *tot_energies;
  double *forces;
  double *forces_single;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  N = 1024;
  int particles_per_process = floor(N/size);
  int start_index = rank * particles_per_process;
  int end_index = start_index + particles_per_process;
  int residual_particles = N - (particles_per_process * size);
  int index_residual = N - (residual_particles + 1);

  r_c = 10.0;
  h = 1e-20;
  positions = (double *)malloc(N * 3 * sizeof(double));
  distances = (double *)malloc(N * N * sizeof(double));
  potentials = (double *)malloc(N * N * sizeof(double));
  tot_energies = (double *)malloc(N * sizeof(double));
  forces = (double *)malloc(N * N * 3 * sizeof(double));
  forces_single = (double *)malloc(N * 3 * sizeof(double));

  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0){    
    read_positions(positions, N);
    calculate_distances(positions, distances, N);

    start_time = MPI_Wtime();
    for (int i = 0; i < N; ++i) {
        for (int direction = 0; direction < 3; ++direction) {
            forces_single[i * 3 + direction] = calc_force_per_particle(positions, distances, N, i, direction);
        }
    }
    end_time = MPI_Wtime();
    double time_spent = end_time - start_time;
    printf("Sequential calculation:\n");
    printf("Time spent for calculating the forces of %d particles is %f seconds. From rank: %d, of %d processors\n", N, time_spent, rank, size);
//    print_forces_single(forces_single, N);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(rank == 0){
    printf("Starting  parallel calculation\n");
    printf("There are %d particles per process and %d residual particles. \n The index of the first resudial particle is %d\n", particles_per_process, residual_particles, index_residual);
    start_time2 = MPI_Wtime();
  }
  for (int i = start_index; i < end_index; ++i) {
      for (int direction = 0; direction < 3; ++direction) {
          forces_single[i * 3 + direction] = calc_force_per_particle(positions, distances, N, i, direction);
      }
  }
  if(rank == 0){
    for (int i = index_residual; i < N; ++i) {
      for (int direction = 0; direction < 3; ++direction) {
          forces_single[i * 3 + direction] = calc_force_per_particle(positions, distances, N, i, direction);
      }
    }
  }
  if(rank == 0){
    end_time2 = MPI_Wtime();
    double time_spent2 = end_time2 - start_time2;
    printf("Parallel calculation:\n");
    printf("Time spent for calculating the forces of %d particles is %f seconds. From rank: %d, of %d processors\n", N, time_spent2, rank, size);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  printf("Hello from %d of %d processors.\n", rank, size);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  free(positions);
  free(distances);
  free(potentials);
  free(tot_energies);
  free(forces);
  free(forces_single);
  return 0;
}

