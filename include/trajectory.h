#pragma once

#include "geometry.h"
#include     <chrono>
#include    <fstream>
#include    <sstream>

class Trajectory {
public:
    static Trajectory Load(const std::string& movie);
    void move(const glm::vec3& vector);
    void render(const Shader& shader);

private:
    std::chrono::high_resolution_clock::time_point timestamp; std::vector<Geometry> geoms; bool paused = false; float wait = 15.997; int frame = 0;
};
