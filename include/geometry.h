#pragma once

#include                "glfwpointer.h"
#include                   "constant.h"
#include                       "mesh.h"
#include                 <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include                      <sstream>

class Geometry {
    struct Object {
        glm::mat4 get_model(glm::mat4 s = glm::mat4(1)) const {
            return translate * rotate * s * scale;
        }
        glm::vec3 get_position() const {
            return glm::vec3(translate[3]);
        }
        glm::mat4 translate, rotate, scale; std::string name;
    };

public:
    static Geometry Load(std::stringstream& file);
    glm::vec3 get_center() const;
    void move(const glm::vec3& vector);
    void render(const Shader& shader) const;
    void rebind(float factor);

    inline static std::unordered_map<std::string, Mesh> meshes;

private:
    std::vector<Object> objects;
};
