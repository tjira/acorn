#include "shader.h"

Shader::Shader(const std::string& vertex, const std::string& fragment) : id(glCreateProgram()) {
    // extract the raw strings
    const char *fsCode = fragment.c_str(), *vsCode = vertex.c_str();

    // create the shaders
    unsigned int fs = glCreateShader(GL_FRAGMENT_SHADER);
    unsigned int vs = glCreateShader(GL_VERTEX_SHADER);

    // set the shader sources
    glShaderSource(vs, 1, &vsCode, nullptr);
    glShaderSource(fs, 1, &fsCode, nullptr);

    // compile the shaders
    glCompileShader(vs), glCompileShader(fs);

    // check for errors
    check_error(vs, "vertex"), check_error(fs, "fragment");

    // attach shaders to the program
    glAttachShader(id, vs), glAttachShader(id, fs);

    // link and validate the program
    glLinkProgram(id), glValidateProgram(id);

    // detach shaders from program
    glDetachShader(id, vs), glDetachShader(id, fs);

    // delete the shaders
    glDeleteShader(vs), glDeleteShader(fs), use();
}

Shader::~Shader() {
    glDeleteProgram(id);
}

void Shader::check_error(unsigned int shader, const std::string& title) const {
    // extract status code
    int success; glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    // extract message if fail
    if(int length = 0; success == 0) {

        // extract message length
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);

        // extract the message
        std::vector<char> log(length); glGetShaderInfoLog(shader, length, &length, &log[0]);

        // throw the error
        throw std::runtime_error(title + ": " + std::string(log.begin(), log.end()));
    }
}

void Shader::use() const {
    glUseProgram(id);
}

template <typename T>
void Shader::set(const std::string& name, T value) const {
    if constexpr (std::is_same<T, int>()      )        glUniform1i(glGetUniformLocation(id, name.c_str()), value                                 );
    if constexpr (std::is_same<T, float>()    )        glUniform1f(glGetUniformLocation(id, name.c_str()), value                                 );
    if constexpr (std::is_same<T, glm::vec3>())        glUniform3f(glGetUniformLocation(id, name.c_str()), value[0], value[1], value[2]          );
    if constexpr (std::is_same<T, glm::vec4>())        glUniform4f(glGetUniformLocation(id, name.c_str()), value[0], value[1], value[2], value[3]);
    if constexpr (std::is_same<T, glm::mat4>()) glUniformMatrix4fv(glGetUniformLocation(id, name.c_str()), 1,        GL_FALSE, &value[0][0]      );
}

template void Shader::set<float>    (const std::string& name, float value    ) const;
template void Shader::set<glm::vec3>(const std::string& name, glm::vec3 value) const;
template void Shader::set<glm::vec4>(const std::string& name, glm::vec4 value) const;
template void Shader::set<glm::mat4>(const std::string& name, glm::mat4 value) const;
template void Shader::set<int>      (const std::string& name, int value      ) const;
