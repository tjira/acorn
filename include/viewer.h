#pragma once

#include                 "glfwpointer.h"
#include                   <glad/glad.h>
#include                  <GLFW/glfw3.h>
#include             <ImGuiFileDialog.h>
#include    <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>
#include                      <implot.h>

class Viewer {
public:
    Viewer(const std::vector<std::string>& inputs);
    void gui() const;

private:
    GLFWPointer pointer; std::vector<std::string> inputs;
};
