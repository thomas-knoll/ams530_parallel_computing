import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualize_particles(file_name, output_file):
    # Read particle coordinates from the file
    with open(file_name, 'r') as file:
        particles = [list(map(float, line.split())) for line in file]

    particles = sorted(particles, key=lambda x: x[0])  # Sort particles by index

    # Extract x, y, and z coordinates
    x = [particle[1] for particle in particles]
    y = [particle[2] for particle in particles]
    z = [particle[3] for particle in particles]

    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Create box
    box_edges = [(-40, -40, -40), (40, -40, -40), (40, -40, 40), (-40, -40, 40), (-40, -40, -40),
                 (-40, 40, -40), (40, 40, -40), (40, 40, 40), (-40, 40, 40), (-40, 40, -40),
                 (-40, -40, -40), (-40, 40, -40), (40, 40, -40), (40, -40, -40), (40, -40, 40),
                 (40, 40, 40), (-40, 40, 40), (-40, -40, 40), (-40, -40, -40)]

    box_edges = list(zip(*box_edges))
    ax.plot(box_edges[0], box_edges[1], box_edges[2], color='black')

    # Visualize particles as darkgreen dots
    ax.scatter(x, y, z, c='darkgreen', marker='.')

    # Add labels at the middle of three edges in each dimension
#    ax.text(0, -40, 0, '80', color='black', fontsize=8)
#    ax.text(-40, 0, 0, '80', color='black', fontsize=8)
#    ax.text(0, 0, -40, '80', color='black', fontsize=8)

    # Set axis limits
    ax.set_xlim(-40, 40)
    ax.set_ylim(-40, 40)
    ax.set_zlim(-40, 40)

    # Hide axes
    ax.set_axis_off()
    
    # Rotate the view
    ax.view_init(elev=15, azim=30)  # Set the elevation and azimuthal angles

    # Save plot as PNG
    plt.savefig(output_file, bbox_inches='tight')
#    plt.show()

if __name__ == "__main__":
    input_file = "particles_1024.txt"
    output_file = "particles_1024.png"

    visualize_particles(input_file, output_file)
    print(f"Particle visualization saved to {output_file}.")

