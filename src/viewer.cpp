#include "viewer.h"

std::string vertex = R"(
#version 420 core

layout(location = 0) in vec3 i_position;
layout(location = 1) in vec3 i_normal;
layout(location = 2) in vec3 i_color;

uniform mat4 u_model, u_view, u_proj;

out mat3 transform; out vec3 fragment, normal, color;

void main() {
    normal = normalize(mat3(transpose(inverse(u_model))) * i_normal);
    fragment = vec3(u_model * vec4(i_position, 1)), color = i_color;
    gl_Position = u_proj * u_view * vec4(fragment, 1);
    transform = inverse(mat3(u_view));
})";

std::string fragment = R"(
#version 420 core

struct Light {vec3 position; float ambient, diffuse, specular, shininess;};

uniform Light u_light; uniform vec3 u_camera;

in mat3 transform; in vec3 fragment, normal, color;

out vec4 o_color;

void main() {
    vec3 lightPos = transform * u_light.position, reflection = reflect(-normalize(lightPos), normal), direction = normalize(u_camera - fragment);
    vec3 specular = vec3(pow(max(dot(direction, reflection), 0), u_light.shininess)),  diffuse = vec3(max(dot(normal, normalize(lightPos)), 0));
    o_color = vec4((vec3(u_light.ambient) + u_light.diffuse * diffuse + u_light.specular * specular), 1) * vec4(color, 1);
})";

void set(const Shader& shader, const GLFWPointer::Camera& camera, const GLFWPointer::Light& light) {
    // use the shader
    shader.use();

    // set the camera
    shader.set<glm::vec3>("u_camera", -glm::inverse(glm::mat3(camera.view)) * glm::vec3(camera.view[3]));

    // set other variables
    shader.set<glm::mat4>("u_view",            camera.view    );
    shader.set<glm::mat4>("u_proj",            camera.proj    );
    shader.set<glm::vec3>("u_light.position",  light.position );
    shader.set<float>    ("u_light.shininess", light.shininess);
    shader.set<float>    ("u_light.specular",  light.specular );
    shader.set<float>    ("u_light.ambient",   light.ambient  );
    shader.set<float>    ("u_light.diffuse",   light.diffuse  );
}

void resizeCallback(GLFWwindow* window, int width, int height) {
    // extract the pointer and check if the window is present after resizing
    if (GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window); width > 0 && height > 0) {

        // change the perspective
        pointer->camera.proj = glm::perspective(glm::radians(45.0f), (float)width / height, 0.01f, 1000.0f);

        // change the window size variables and the viewport
        pointer->width = width, pointer->height = height; glViewport(0, 0, width, height);
    }
}

void keyCallback(GLFWwindow* window, int key, int, int action, int mods) {
    // extract the pointer and check for button press
    if (GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window); action == GLFW_PRESS) {

        // flag variables
        if      (key == GLFW_KEY_F1   ) pointer->flags.options = !pointer->flags.options;
        else if (key == GLFW_KEY_F2   ) pointer->flags.plot    =    !pointer->flags.plot;
        else if (key == GLFW_KEY_F12  ) pointer->flags.info    =    !pointer->flags.info;
        else if (key == GLFW_KEY_SPACE) pointer->flags.pause   =   !pointer->flags.pause;
        else if (key == GLFW_KEY_F11  ) {

            // define the necessary variables for fullscreen
            static int xpos0, ypos0, width0, height0; int xpos, ypos, width, height;

            // if the monitor is requested to set fullscreen
            if (pointer->flags.fullscreen = !pointer->flags.fullscreen; pointer->flags.fullscreen) {

                // save the initial variables
                glfwGetWindowSize(pointer->window, &width0, &height0), glfwGetWindowPos(pointer->window, &xpos0, &ypos0);

                // get the variables for the active monitor
                glfwGetMonitorWorkarea(glfwGetPrimaryMonitor(), &xpos, &ypos, &width, &height);

                // set the fullscreen mode
                glfwSetWindowMonitor(pointer->window, glfwGetPrimaryMonitor() , 0, 0, width, height, 1.0 / 60);

            // revert to the original state if monitor is already in the fullscreen mode
            } else glfwSetWindowMonitor(pointer->window, nullptr , xpos0, ypos0, width0, height0, 1.0 / 60), resizeCallback(pointer->window, width0, height0);
        }
    }
}

void positionCallback(GLFWwindow* window, double x, double y) {
    // get the glfw pointer
    GLFWPointer* pointer = (GLFWPointer*)glfwGetWindowUserPointer(window);

    // check if the left mouse button is pressed
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS && !ImGui::GetIO().WantCaptureMouse) {

        // calculate the rotation axes
        glm::vec3 xaxis = glm::inverse(glm::mat3(pointer->camera.view)) * glm::vec3(0, 1, 0);
        glm::vec3 yaxis = glm::inverse(glm::mat3(pointer->camera.view)) * glm::vec3(1, 0, 0);

        // rotate the camera
        pointer->camera.view = glm::rotate(pointer->camera.view, 0.01f * ((float)y - pointer->mouse.y), yaxis);
        pointer->camera.view = glm::rotate(pointer->camera.view, 0.01f * ((float)x - pointer->mouse.x), xaxis);
    }
    
    // set the mopuse position
    pointer->mouse = {x, y};
}

void scrollCallback(GLFWwindow* window, double, double dy) {
    if (!ImGui::GetIO().WantCaptureMouse) ((GLFWPointer*)glfwGetWindowUserPointer(window))->camera.view *= glm::mat4(glm::mat3(1.0f + 0.08f * (float)dy));
}

Viewer::Viewer(const std::vector<std::string>& inputs) : inputs(inputs) {
    // initialize GLFW and throw error if failed
    if(!glfwInit()) throw std::runtime_error("ERROR DURING GLFW INITIALIZATION");

    // pass OpenGL version and other hints
    glfwWindowHint(GLFW_OPENGL_PROFILE,        GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, pointer.major           );
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, pointer.minor           );
    glfwWindowHint(GLFW_SAMPLES,               pointer.samples         );

    // create the window
    if (pointer.window = glfwCreateWindow(pointer.width, pointer.height, pointer.title.c_str(), nullptr, nullptr); !pointer.window) throw std::runtime_error("ERROR DURING WINDOW CREATION");

    // initialize GLAD
    if (glfwMakeContextCurrent(pointer.window); !gladLoadGL()) throw std::runtime_error("ERROR DURING GLAD INITIALIZATION");

    // set some GLAD options
    glEnable(GL_DEPTH_TEST), glEnable(GL_CULL_FACE);

    // set some GLFW options
    glfwSetWindowUserPointer(pointer.window, &pointer); glfwSwapInterval(1);

    // set event callbacks
    glfwSetCursorPosCallback(pointer.window, positionCallback);
    glfwSetWindowSizeCallback(pointer.window, resizeCallback);
    glfwSetScrollCallback(pointer.window, scrollCallback);
    glfwSetKeyCallback(pointer.window, keyCallback);

    // initialize GUI
    ImGui::CreateContext(); ImPlot::CreateContext(); ImGui_ImplOpenGL3_Init("#version 420"); ImGui_ImplGlfw_InitForOpenGL(pointer.window, true); ImGui::GetIO().IniFilename = nullptr;

    // initialize camera matrices
    pointer.camera.proj = glm::perspective(glm::radians(45.0f), (float)pointer.width / pointer.height, 0.01f,                       1000.0f);
    pointer.camera.view = glm::lookAt     ({0.0f, 0.0f, 5.0f},  glm::vec3(0.0f),                       glm::vec3(0.0f, 1.0f, 0.0f)         );

    {
        // read the input matrices
        std::vector<Eigen::MatrixXd> matrices; for (const std::string& input : inputs) if (input.substr(input.length() - 3) == "mat") matrices.push_back(read(input));

        // read the input molecules
        std::vector<Trajectory> trajectories; for (const std::string& input : inputs) if (input.substr(input.length() - 3) == "xyz") trajectories.push_back(Trajectory::Load(input));

        // set the default indices
        std::vector<std::tuple<int, std::vector<int>, int>> indices; for (size_t i = 0; i < matrices.size(); i++) indices.push_back({1, {1}, 0});

        // create the meshes of atoms
        for (auto& [symbol, object] : ptable) Geometry::meshes[symbol] = Mesh::Icosphere(SUBDIVISIONS, SMOOTH, symbol), Geometry::meshes.at(symbol).setColor(object.color);
        
        // create the cyllinder mesh for bond
        Geometry::meshes["bond"] = Mesh::Cylinder(SECTORS, SMOOTH, "bond");

        // create the shaders
        Shader shader(vertex, fragment);

        // enter the render loop
        while (!glfwWindowShouldClose(pointer.window)) {

            // clear the color and depth buffer
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // set shader variables
            set(shader, pointer.camera, pointer.light);

            // render stuff
            for (Trajectory& traj : trajectories) {traj.render(shader);} gui(matrices, indices); 
            
            // swap buffers and poll events
            glfwSwapBuffers(pointer.window); glfwPollEvents();
        }
    }

    // destroy the GUI
    ImGui_ImplGlfw_Shutdown(); ImGui_ImplOpenGL3_Shutdown(); ImPlot::DestroyContext(); ImGui::DestroyContext();
}

void Viewer::gui(const std::vector<Eigen::MatrixXd>& matrices, std::vector<std::tuple<int, std::vector<int>, int>>& indices) {
    // begin the frame and define id
    ImGui_ImplOpenGL3_NewFrame(); ImGui_ImplGlfw_NewFrame(); ImGui::NewFrame(); int id = 0;

    // plotter window
    if (pointer.flags.plot) {

        // begin the plotter window
        ImGui::SetNextWindowPos({0, 0}, ImGuiCond_Once); ImGui::Begin("Plotter", &pointer.flags.plot);

        // begin the plot
        if (ImPlot::BeginPlot("##plot")) {

            // plot the data
            for (size_t i = 0; i < matrices.size(); i++) for (size_t j = 0; j < std::get<1>(indices.at(i)).size(); j++) {
                ImPlot::PlotLine(("##" + std::to_string(id++)).c_str(), matrices.at(i).col(0).data(), matrices.at(i).col(std::get<1>(indices.at(i)).at(j)).data(), matrices.at(i).rows());
            }

            // end the plot
            ImPlot::EndPlot();
        }

        // show the matrix plot options
        for (size_t i = 0; i < indices.size(); i++) {

            // add the column sliders
            for (size_t j = 0; j < std::get<1>(indices.at(i)).size(); j++) {
                ImGui::SetNextItemWidth(50); ImGui::PushID(id++); if (j) {ImGui::SameLine();} ImGui::SliderInt("##", &std::get<1>(indices.at(i)).at(j), 0, matrices.at(i).cols() - 1); ImGui::PopID();
            }

            // add the button to add a column to the plot
            ImGui::SameLine(); if (ImGui::Button("Add Column")) std::get<1>(indices.at(i)).push_back(std::get<1>(indices.at(i)).at(std::get<1>(indices.at(i)).size() - 1) + 1);

            // add the button to remove a column from the plot
            ImGui::SameLine(); if (ImGui::Button("Remove Column")) if (std::get<1>(indices.at(i)).size() > 1) std::get<1>(indices.at(i)).pop_back();

            // add the animation step slider
            ImGui::SameLine(); ImGui::SetNextItemWidth(50); ImGui::SliderInt(("##" + std::to_string(id++)).c_str(), &std::get<2>(indices.at(i)), 0, matrices.at(i).cols());
            
            // increase the indices
            for (int& index : std::get<1>(indices.at(i))) if (!pointer.flags.pause) index = (index + std::get<2>(indices.at(i))) % (matrices.at(i).cols() - std::get<0>(indices.at(i)));
        }

        // end the plotter window
        ImGui::End();
    }

    // end the frame
    ImGui::Render(); ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

Eigen::MatrixXd Viewer::read(const std::string& path) const {
    // open the file and extract shape
    std::ifstream fstream(path); int rows, cols; fstream >> rows >> cols;

    // create the matrix and fill it
    Eigen::MatrixXd matrix(rows, cols); for (int i = 0; i < rows; i++) for (int j = 0; j < cols; j++) fstream >> matrix(i, j);

    // return the matrix
    return matrix;
}
