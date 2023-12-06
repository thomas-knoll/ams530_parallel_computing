import random

def generate_particles(num_particles, file_name):
    with open(file_name, 'w') as file:
        for i in range(num_particles):
            x = random.uniform(-40, 40)
            y = random.uniform(-40, 40)
            z = random.uniform(-40, 40)
            file.write(f"{i} {x} {y} {z}\n")

if __name__ == "__main__":
    num_particles = 2**10  # Number of particles (N = 1024)
#    num_particles = 300  # Number of particles (N = 1024)
    output_file = "particles_1024.txt"

    generate_particles(num_particles, output_file)
    print(f"Particles generated and saved to {output_file}.")

