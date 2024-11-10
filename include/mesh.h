#pragma once

#include  "buffer.h"
#include  "shader.h"
#include <algorithm>

class Mesh {
public:
    Mesh() {}; Mesh(std::vector<Vertex> data, const std::string& name = "mesh") : name(name), model(1.0f), buffer(data) {};
    static Mesh Cylinder(int sectors, bool smooth, const std::string& name = "cylinder");
    static Mesh Icosphere(int subdivisions, bool smooth, const std::string& name = "icosphere");
    std::string getName() const;
    glm::vec3 getPosition() const;
    void render(const Shader& shader, const glm::mat4& transform = glm::mat4(1.0f)) const;
    void setColor(const glm::vec3& color);
    void setModel(const glm::mat4& model);

private:
    std::string name; glm::mat4 model; Buffer buffer;
};
