#include "trajectory.h"

Trajectory Trajectory::Load(const std::string& filename) {
    // define the trajectory and line container and open the file and skip first line
    Trajectory trajectory; std::ifstream file(filename); std::string line; std::getline(file, line);

    // extract the length of the geometry and number of geometries
    int length = std::stoi(line), block = 1; while (std::getline(file, line)) block++;

    // create the vector of geometris and reset the file reader.
    trajectory.geoms.resize(block / (length + 2)); file.clear(), file.seekg(0);

    // read the individual geometries.
    for (int i = 0; i < block / (length + 2); i++) {
        std::stringstream ss; for (int j = 0; j < length + 2; j++) {std::getline(file, line); ss << (j ? "\n" : "" ) << line;} trajectory.geoms.at(i) = Geometry::Load(ss);
    }

    // set the initialization timestamp, center the trajectory and return the trajectory
    trajectory.timestamp = std::chrono::high_resolution_clock().now(), trajectory.move(-trajectory.geoms.at(0).get_center()); return trajectory;
}

void Trajectory::move(const glm::vec3& vector) {
    for (Geometry& geom : geoms) geom.move(vector);
}

void Trajectory::render(const Shader& shader) {
    // if geoms not empty
    if (geoms.size()) {

        // measure the elapsed time from the last render
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock().now() - timestamp).count();

        // if the elapsed time is longer than the waiting period, increase the frame index
        if (elapsed > wait && !paused && wait > 0) frame = (frame + (int)(elapsed / wait)) % (int)geoms.size(), timestamp = std::chrono::high_resolution_clock().now();

        // render the geometry
        geoms.at(frame).render(shader);
    }
}
