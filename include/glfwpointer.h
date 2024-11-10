#pragma once

#include <glm.hpp>
#include  <string>

#define WIDTH 1024
#define HEIGHT 576
#define SUBDIVISIONS 2
#define SECTORS 16
#define SMOOTH 1
#define BINDINGFACTOR 0.013
#define BONDSIZE 0.09
#define ATOMSIZEFACTOR 0.007

struct GLFWwindow;

struct GLFWPointer {
    std::string title = "Quantum Acorn"; glm::vec2 mouse; GLFWwindow* window;
    int width = WIDTH, height = HEIGHT, samples = 16, major = 4, minor = 2;
    struct Camera {
        glm::mat4 view, proj;
    } camera;
    struct Flags {
        bool fullscreen = false, info = false, options = false, pause = false, plot = false, system = false, ptable = false;
    } flags;
    struct Light {
        float ambient = 0.4f, diffuse = 0.2f, specular = 0.4f, shininess = 4.0f;
        glm::vec3 position = {1.0f, 1.0f, 1.0f};
    } light;
};
