#include "viewer.h"

#include <iostream>

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

    // initialize GUI
    ImGui::CreateContext(); ImPlot::CreateContext(); ImGui_ImplOpenGL3_Init("#version 420"); ImGui_ImplGlfw_InitForOpenGL(pointer.window, true); ImGui::GetIO().IniFilename = nullptr;

    {
        // read the inputs
        std::vector<Eigen::MatrixXd> matrices; for (const std::string& input : inputs) matrices.push_back(read(input));

        // set the default indices
        std::vector<std::tuple<int, std::vector<int>, int>> indices; for (int i = 0; i < matrices.size(); i++) indices.push_back({1, {1}, 0});

        // enter the render loop
        while (!glfwWindowShouldClose(pointer.window)) {

            // clear the color and depth buffer
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // render stuff
            gui_plot(matrices, indices);
            
            // swap buffers and poll events
            glfwSwapBuffers(pointer.window); glfwPollEvents();
        }
    }

    // destroy the GUI
    ImGui_ImplGlfw_Shutdown(); ImGui_ImplOpenGL3_Shutdown(); ImPlot::DestroyContext(); ImGui::DestroyContext();
}

void Viewer::gui_plot(const std::vector<Eigen::MatrixXd>& matrices, std::vector<std::tuple<int, std::vector<int>, int>>& indices) const {
    // begin the frame and extract window size
    ImGui_ImplOpenGL3_NewFrame(); ImGui_ImplGlfw_NewFrame(); ImGui::NewFrame(); int width, height; glfwGetWindowSize(pointer.window, &width, &height);

    // show demo
    // ImGui::ShowDemoWindow(); ImPlot::ShowDemoWindow();

    // begin the plotter window
    ImGui::SetNextWindowPos({0, 0}); ImGui:: SetNextWindowSize({(float)width, (float)height}); ImGui::Begin("Plotter", nullptr);

    // begin the plot
    if (ImPlot::BeginPlot("##plot")) {


        // plot the data
        for (int i = 0; i < matrices.size(); i++) for (int j = 0; j < std::get<1>(indices.at(i)).size(); j++) {
            Eigen::VectorXd x = matrices.at(i).col(0), y = matrices.at(i).col(std::get<1>(indices.at(i)).at(j)); ImPlot::PlotLine("##", x.data(), y.data(), x.rows());
        }

        // end the plot
        ImPlot::EndPlot();
    }

    // show the matrix plot options
    for (size_t i = 0, j = 0; i < indices.size(); i++) {

        // add the column sliders
        for (size_t k = 0; k < std::get<1>(indices.at(i)).size(); k++, j++) {
            ImGui::SetNextItemWidth(100); ImGui::PushID(j); if (k) {ImGui::SameLine();} ImGui::SliderInt("##", &std::get<1>(indices.at(i)).at(k), 0, matrices.at(i).cols() - 1); ImGui::PopID();
        }

        // add the button to add a column to the plot
        ImGui::SameLine(); if (ImGui::Button("Add Column")) std::get<1>(indices.at(i)).push_back(std::get<1>(indices.at(i)).at(std::get<1>(indices.at(i)).size() - 1) + 1);

        // add the button to remove a column from the plot
        ImGui::SameLine(); if (ImGui::Button(("Remove Column"))) if (std::get<1>(indices.at(i)).size() > 1) std::get<1>(indices.at(i)).pop_back();

        // add the animation step slider
        ImGui::SameLine(); ImGui::SetNextItemWidth(100); ImGui::PushID(j); ImGui::SliderInt(("##" + std::to_string(i)).c_str(), &std::get<2>(indices.at(i)), 0, matrices.at(i).cols() - 1); ImGui::PopID();
        
        // increase the indices
        for (int& index : std::get<1>(indices.at(i))) index = (index + std::get<2>(indices.at(i))) % matrices.at(i).cols();
    }

    // end the plotter window
    ImGui::End();

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
