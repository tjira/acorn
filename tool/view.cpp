#include <glad/gl.h>
#include <GLFW/glfw3.h>
#include <fstream>
#include <argparse/argparse.hpp>
#include <glm/gtc/matrix_transform.hpp>

struct GLFWPointer {
    std::string title="Acorn"; glm::vec2 mouse; GLFWwindow* window;
    int width=1024, height=576, samples=16, major=4, minor=2;
    struct Camera {
        glm::mat4 view, proj;
    } camera{};
    struct Light {
        glm::vec3 position={1.0f, 1.0f, 1.0f}; float ambient=0.4f, diffuse=0.2f, specular=0.4f, shininess=4.0f;
    } light{};
};

struct Vertex {
    glm::vec3 position, normal=glm::vec3(0), color=glm::vec3(1);
};

class Buffer {
public:
    // buffer constructor
    Buffer(const std::vector<Vertex>& data);

    // getters
    std::vector<Vertex> getData() const {return data;}
    int getSize() const {return data.size();}

    // state functions
    void bind() const {glBindVertexArray(vao);}

private:
    std::vector<Vertex> data;
    unsigned int vao, vbo;
};

class Shader {
public:
    // constructor and destructor
    ~Shader() {glDeleteProgram(id);} Shader(const std::string& vertex, const std::string& fragment);

    // setters
    template <typename T> void set(const std::string& name, T value) const;

    // state functions
    void use() const {glUseProgram(id);}

private:
    unsigned int id;
};

class Mesh {
public:
    // mesh constructor
    Mesh(std::vector<Vertex> data) : model(1.0f), buffer(data) {};

    // static constructors
    static std::vector<std::vector<Mesh>> Wavefunctions(const std::string& path, int dim);
    static std::vector<Mesh> Matrix(const std::string& path, int dim);

    // setters
    void setColor(const glm::vec3& color); void setModel(const glm::mat4& model);

    // state functions
    void render(const Shader& shader, const glm::mat4& transform=glm::mat4(1.0f)) const;

private:
    Buffer buffer; glm::mat4 model;
};

Buffer::Buffer(const std::vector<Vertex>& data) : data(data) {
    glGenVertexArrays(1, &vao), glGenBuffers(1, &vbo), glBindBuffer(GL_ARRAY_BUFFER, vbo), glBindVertexArray(vao);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, color));
    glEnableVertexAttribArray(0), glEnableVertexAttribArray(1), glEnableVertexAttribArray(2);
    glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(Vertex), data.data(), GL_STATIC_DRAW);
}

std::vector<Mesh> Mesh::Matrix(const std::string& path, int dim) {
    // define the data and the file stream
    std::vector<std::vector<Vertex>> data; std::ifstream file(path); std::string line;

    while (std::getline(file, line)) {
        // define the line stream, variables and index
        std::stringstream stream(line); std::vector<float> xyz(3, 0.0f); int i = 0;

        // fill the independent variables
        for (int j = 0; j < dim; j++) stream >> xyz.at(j);

        while (stream >> xyz.at(dim)) {
            // push the mesh data
            if (data.size() <= i) {data.push_back({});}

            // push the current point and increase index
            for (int j = 0; j < 2; j++) {data.at(i).push_back({{xyz.at(0), xyz.at(1), xyz.at(2)}});} i++;
        }
    }

    // remove one first and one last element
    for (auto& set : data) set.erase(set.begin()), set.erase(set.end() - 1);

    // create the mesh vector and return
    std::vector<Mesh> mesh; for (const auto& set : data) {mesh.emplace_back(set);} return mesh;
}

std::vector<std::vector<Mesh>> Mesh::Wavefunctions(const std::string& path, int dim) {
    // define the data, file stream and some containers
    std::vector<std::vector<std::vector<Vertex>>> data; std::vector<std::vector<Mesh>> mesh;
    std::ifstream file(path); std::string line; std::vector<float> r;

    // read the header
    std::getline(file, line);

    // read the independent variables
    for (int i = 0; i < dim; i++) {
        std::getline(file, line); std::stringstream stream(line);
        float value; while (stream >> value) r.push_back(value);
    }

    while (std::getline(file, line)) {
        // push the time data
        data.push_back({});
        
        // for every wfn data
        for (int i = 0; i < 2 * dim; i++) {
            // push the wfn data and get the line
            data.back().push_back({}); std::getline(file, line);

            // create the stream and the point container
            std::stringstream stream(line); float value;

            // fill the data
            while (stream >> value) {
                if (dim == 2) data.back().at(i).push_back({{r.at(data.back().at(i).size() % r.size()), r.at(data.back().at(i).size() / r.size()), value}});
                else if (dim == 1) data.back().at(i).push_back({{r.at(data.back().at(i).size() % r.size()), value, 0}});
            }
        }
    }

    // fill the mesh vector and return
    for (const auto& time : data) {std::vector<Mesh> meshes;
        for (const auto& set : time) meshes.emplace_back(set); mesh.push_back(meshes);
    } return mesh;
}

void Mesh::render(const Shader& shader, const glm::mat4& transform) const {
    shader.use(), shader.set<glm::mat4>("u_model", transform * model);
    buffer.bind(), glDrawArrays(GL_POINTS, 0, buffer.getSize());
}

Shader::Shader(const std::string& vertex, const std::string& fragment) : id(glCreateProgram()) {
    // create the shader pointers
    const char *fsCode = fragment.c_str(), *vsCode = vertex.c_str();

    // create the shaders
    unsigned int fs = glCreateShader(GL_FRAGMENT_SHADER);
    unsigned int vs = glCreateShader(GL_VERTEX_SHADER);

    // set the shader sources
    glShaderSource(vs, 1, &vsCode, nullptr);
    glShaderSource(fs, 1, &fsCode, nullptr);

    // compile the shader
    glCompileShader(vs), glCompileShader(fs);

    // attach and link the shader
    glAttachShader(id, vs), glAttachShader(id, fs);
    glLinkProgram(id), glValidateProgram(id);

    // detach, delete and use the shader
    glDetachShader(id, vs), glDetachShader(id, fs);
    glDeleteShader(vs), glDeleteShader(fs), use();
}

template <typename T>
void Shader::set(const std::string& name, T value) const {
    if constexpr (std::is_same<T, int>()) glUniform1i(glGetUniformLocation(id, name.c_str()), value);
    if constexpr (std::is_same<T, float>()) glUniform1f(glGetUniformLocation(id, name.c_str()), value);
    if constexpr (std::is_same<T, glm::vec3>()) glUniform3f(glGetUniformLocation(id, name.c_str()), value[0], value[1], value[2]);
    if constexpr (std::is_same<T, glm::vec4>()) glUniform4f(glGetUniformLocation(id, name.c_str()), value[0], value[1], value[2], value[3]);
    if constexpr (std::is_same<T, glm::mat4>()) glUniformMatrix4fv(glGetUniformLocation(id, name.c_str()), 1, GL_FALSE, &value[0][0]);
}

void keyCallback(GLFWwindow* window, int key, int, int action, int mods) {
    if (GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window); action == GLFW_PRESS) {}
}

void positionCallback(GLFWwindow* window, double x, double y) {
    // get the pointer
    GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window);

    // check if the left mouse button is pressed
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
        // get the x and y axis of rotation
        glm::vec3 xaxis = glm::inverse(glm::mat3(pointer->camera.view)) * glm::vec3(1, 0, 0);
        glm::vec3 yaxis = glm::inverse(glm::mat3(pointer->camera.view)) * glm::vec3(0, 1, 0);
        
        // update the view matrix
        pointer->camera.view = glm::rotate(pointer->camera.view, 0.01f * ((float)y - pointer->mouse.y), xaxis);
        pointer->camera.view = glm::rotate(pointer->camera.view, 0.01f * ((float)x - pointer->mouse.x), yaxis);
    } else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
        // get the x and y axis of traanslation
        glm::vec3 xaxis = glm::inverse(glm::mat3(pointer->camera.view)) * glm::vec3(1, 0, 0);
        glm::vec3 yaxis = glm::inverse(glm::mat3(pointer->camera.view)) * glm::vec3(0, 1, 0);
        
        // update the view matrix
        pointer->camera.view = glm::translate(pointer->camera.view, 0.01f * ((float)x - pointer->mouse.x) * xaxis);
        pointer->camera.view = glm::translate(pointer->camera.view, 0.01f * (pointer->mouse.y - (float)y) * yaxis);
    }

    // update the mouse position
    pointer->mouse = {x, y};
}

void resizeCallback(GLFWwindow* window, int width, int height) {
    // get the pointer and check if the width and height are valid
    if (GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window); width > 0 && height > 0) {
        // update the projection matrix
        pointer->camera.proj = glm::perspective(glm::radians(45.0f), (float)width / height, 1e-2f, 1e3f);

        // update the viewport and the window size
        pointer->width = width, pointer->height = height; glViewport(0, 0, width, height);
    }
}

void scrollCallback(GLFWwindow* window, double, double dy) {
    // get the pointer and update the view matrix
    ((GLFWPointer*)glfwGetWindowUserPointer(window))->camera.view *= glm::mat4(glm::mat3(1.0f + 0.08f * (float)dy));
}

std::string vertex = R"(
#version 420 core

// define the variables of the passed data
layout(location = 0) in vec3 i_position; layout(location = 1) in vec3 i_normal; layout(location = 2) in vec3 i_color;

// define the set uniforms
uniform mat4 u_model, u_view, u_proj;

// define the output variables
out vec3 fragment, normal, color; out mat3 transform;

void main() {
    // calculate the normal vector
    normal = normalize(mat3(transpose(inverse(u_model))) * i_normal);

    // calculate the fragment position
    fragment = vec3(u_model * vec4(i_position, 1)), color = i_color;

    // calculate the position with the projection and view matrix
    gl_Position = u_proj * u_view * vec4(fragment, 1);

    // assign the view matrix transform
    transform = inverse(mat3(u_view));
})";

std::string fragment = R"(
#version 420 core

// define the light struct
struct Light {vec3 position; float ambient, diffuse, specular, shininess;};

// define the set uniforms
uniform Light u_light; uniform vec3 u_camera;

// define the input variables
in vec3 fragment, normal, color; in mat3 transform;

// output color
out vec4 o_color;

void main() {
    // calculate the light position, reflection vector and  direction vector
    vec3 lightPos = transform * u_light.position, reflection = reflect(-normalize(lightPos), normal), direction = normalize(u_camera - fragment);

    // calculate the specular and diffuse light
    vec3 specular = vec3(pow(max(dot(direction, reflection), 0), u_light.shininess)), diffuse = vec3(max(dot(normal, normalize(lightPos)), 0));

    // calculate the output color from the contributions and add the ambient light
    o_color = vec4((vec3(u_light.ambient) + u_light.diffuse * diffuse + u_light.specular * specular), 1) * vec4(color, 1);
})";

int main(int argc, char** argv) {
    // initialize the argument parser and container for the arguments
    argparse::ArgumentParser program("Acorn View", "1.0", argparse::default_arguments::none);

    // add options to the parser
    program.add_argument("input").help("Input file.").nargs(argparse::nargs_pattern::any);
    program.add_argument("-h").help("-- Display this help message and exit.").default_value(false).implicit_value(true);
    program.add_argument("-d", "--dimension").help("-- Dimensionality of the input.").default_value(1).scan<'i', int>();
    program.add_argument("-s", "--scale").help("-- Scale the input.").default_value(1.0).scan<'g', double>();

    // extract the variables from the command line
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error &error) {
        std::cerr << error.what() << std::endl << std::endl << program; return EXIT_FAILURE;
    }

    // print help if the help flag was provided
    if (program.get<bool>("-h")) {
        std::cout << program.help().str(); return EXIT_SUCCESS;
    }

    // initialize GLFW and throw error if failed
    if(!glfwInit()) {
        throw std::runtime_error("ERROR DURING GLFW INITIALIZATION.");
    }

    // create GLFW variable struct
    GLFWPointer pointer; 

    // pass OpenGL version and other hints
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, pointer.major);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, pointer.minor);
    glfwWindowHint(GLFW_SAMPLES, pointer.samples);

    // create the window
    if (pointer.window = glfwCreateWindow(pointer.width, pointer.height, pointer.title.c_str(), nullptr, nullptr); !pointer.window) {
        throw std::runtime_error("ERROR DURING WINDOW CREATION.");
    }

    // initialize GLAD
    if (glfwMakeContextCurrent(pointer.window); !gladLoadGL(glfwGetProcAddress)) {
        throw std::runtime_error("ERROR DURING GLAD INITIALIZATION.");
    }

    // enable some options
    glfwSetWindowUserPointer(pointer.window, &pointer);
    glfwSwapInterval(1); glEnable(GL_DEPTH_TEST);

    // Set event callbacks
    glfwSetCursorPosCallback(pointer.window, positionCallback);
    glfwSetWindowSizeCallback(pointer.window, resizeCallback);
    glfwSetScrollCallback(pointer.window, scrollCallback);
    glfwSetKeyCallback(pointer.window, keyCallback);

    // 1D line width and point size
    glLineWidth(4), glPointSize(4);

    // initialize camera matrices
    pointer.camera.proj = glm::perspective(glm::radians(45.0f), (float)pointer.width / pointer.height, 0.01f, 1000.0f);
    pointer.camera.view = glm::lookAt({0.0f, 0.0f, 20.0f}, glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));

    {
        // define the containers for the meshes
        std::vector<std::vector<Mesh>> matrices; std::vector<std::vector<std::vector<Mesh>>> wfns;

        // load provided matrices and wfns
        for (const std::string& path : program.get<std::vector<std::string>>("input")) {
            if(path.substr(path.find_last_of(".") + 1) == "mat") {
                matrices.push_back(Mesh::Matrix(path, program.get<int>("-d")));
            } else if(path.substr(path.find_last_of(".") + 1) == "dat") {
                wfns.push_back(Mesh::Wavefunctions(path, program.get<int>("-d")));
            }
        }

        // compile the shader and define frame index
        Shader shader(vertex, fragment); int i = 0;

        // set the light uniforms
        shader.set<glm::vec3>("u_light.position", pointer.light.position);
        shader.set<float>("u_light.shininess", pointer.light.shininess);
        shader.set<float>("u_light.specular", pointer.light.specular);
        shader.set<float>("u_light.ambient", pointer.light.ambient);
        shader.set<float>("u_light.diffuse", pointer.light.diffuse);

        while (!glfwWindowShouldClose(pointer.window)) {
            // clear the color and depth buffer
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // set the camera uniforms
            shader.set<glm::vec3>("u_camera", -glm::inverse(glm::mat3(pointer.camera.view)) * glm::vec3(pointer.camera.view[3]));
            shader.set<glm::mat4>("u_view", pointer.camera.view); shader.set<glm::mat4>("u_proj", pointer.camera.proj);

            // render the meshes
            if (wfns.size()) for (const std::vector<std::vector<Mesh>>& file : wfns) for (const Mesh& mesh : file.at(i)) mesh.render(shader, glm::scale(glm::mat4(1.0f), {1, program.get<double>("-s"), 1}));
            if (matrices.size()) for (const std::vector<Mesh>& file : matrices) for (const Mesh& mesh : file) mesh.render(shader, glm::scale(glm::mat4(1.0f), {1, program.get<double>("-s"), 1}));
            
            // swap buffers, poll events and increase the frame index
            glfwSwapBuffers(pointer.window), glfwPollEvents(); if (wfns.size()) i = (i + 1) % wfns.at(0).size();
        }
    }

    // terminate GLFW
    glfwTerminate();
}
