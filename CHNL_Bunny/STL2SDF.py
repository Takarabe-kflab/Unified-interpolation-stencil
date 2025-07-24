#STL2SDF.py
#Last update: 2025/04/08 by Itsuki Takarabe
import numpy as np
import trimesh
from tqdm import tqdm
import csv

def compute_sdf(mesh, points, chunk_size=10000):
    # compute the distance from the surface to points
    # negative value if the point is inside
    sdf = np.empty(len(points))
    for i in tqdm(range(0, len(points), chunk_size), desc="Calculating SDF"):
        pts_chunk = points[i:i+chunk_size]
        nearest, distance, _ = mesh.nearest.on_surface(pts_chunk)
        inside = mesh.contains(pts_chunk)
        signed_distance = distance.copy()
        signed_distance[inside] = -signed_distance[inside]
        sdf[i:i+chunk_size] = signed_distance
    return sdf

def main():
    # set grid parameters
    lx = 2
    ly = 2
    lz = 2
    nx = 64
    ny = 64
    nz = 64

    # Load STL file
    stl_filename = "stanfordbunny.stl"
    mesh = trimesh.load(stl_filename)
    
    # Put STL file in the center of calculation area
    grid_center = np.array([lx/2, ly/2, lz/2])
    mesh.apply_translation(-mesh.bounding_box.centroid)
    mesh.apply_translation(grid_center)

    # Generate grid coordinate
    x = np.arange(0, lx, lx/nx)
    y = np.arange(0, ly, ly/ny)
    z = np.arange(0, lz, lz/nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T

    # Grid spacing
    x_spacing = lx/nx
    y_spacing = ly/ny
    z_spacing = lz/nz
    # Shift parameters (from the default)
    shifts = [
        ("xshift", np.array([-0.5 * x_spacing, 0, 0])),
        ("yshift", np.array([0, -0.5 * y_spacing, 0])),
        ("zshift", np.array([0, 0, -0.5 * z_spacing])),
        ("noshift", np.array([0, 0, 0]))
    ]

    mesh_centered = mesh.copy()

    # Calculation and output
    for shift_name, offset in shifts:
        mesh_shifted = mesh_centered.copy()
        mesh_shifted.apply_translation(offset)
        print(f"direction: {shift_name}, offset: {offset}")

        sdf = compute_sdf(mesh_shifted, grid_points, chunk_size=10000)
        sdf_grid = sdf.reshape((nx, ny, nz))

        csv_filename = f"sdf/SDF_{shift_name}.csv"
        with open(csv_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # CSV include x,y,z coordinate(based on grid spacing) and SDF
            writer.writerow(["grid_x", "grid_y", "grid_z", "sdf"])
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        writer.writerow([i+1, j+1, k+1, sdf_grid[i, j, k]])
        print(f"output of '{csv_filename}' completed")

if __name__ == "__main__":
    main()
