#pragma once

#include   <glad/glad.h>
#include   <glm/glm.hpp>
#include        <vector>

struct Vertex {
    glm::vec3 position, normal = glm::vec3(0), color = glm::vec3(1);
};

class Buffer {
public:
    Buffer(const Buffer& buffer) : data(buffer.getData()) {generate();};
    Buffer(const std::vector<Vertex>& data) : data(data) {generate();};
    Buffer() : data(0) {generate();}; ~Buffer();
    Buffer& operator=(const Buffer& buffer);
    void bind() const;
    void generate();
    std::vector<Vertex> getData() const { return data; }
    size_t getSize() const { return data.size(); };

private:
    std::vector<Vertex> data;
    unsigned int vao, vbo;
};
