#pragma once

#include <glad/glad.h>
#include <glm/glm.hpp>
#include   <stdexcept>
#include      <vector>

class Shader {
public:
    ~Shader(); Shader(const std::string& vertex, const std::string& fragment);
    void check_error(unsigned int shader, const std::string& title) const;
    template <typename T> void set(const std::string& name, T value) const;
    void use() const;

private:
    unsigned int id;
};
