#pragma once

#include                 "glfwpointer.h"
#include                   <glad/glad.h>
#include                  <GLFW/glfw3.h>
#include             <ImGuiFileDialog.h>
#include    <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>
#include                    <Eigen/Core>
#include                      <implot.h>

class Viewer {
public:
    Viewer(const std::vector<std::string>& inputs);
    void gui_plot(const std::vector<Eigen::MatrixXd>& matrices, std::vector<std::tuple<int, std::vector<int>, int>>& indices) const;
    Eigen::MatrixXd read(const std::string& path) const;

private:
    GLFWPointer pointer; std::vector<std::string> inputs;
};
