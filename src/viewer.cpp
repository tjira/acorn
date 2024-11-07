#include "viewer.h"

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
    glEnable(GL_DEPTH_TEST), glEnable(GL_CULL_FACE), glEnable(GL_STENCIL_TEST);

    // set some GLFW options
    glfwSetWindowUserPointer(pointer.window, &pointer); glfwSwapInterval(1);

    // initialize GUI
    ImGui::CreateContext(); ImPlot::CreateContext(); ImGui_ImplOpenGL3_Init("#version 420"); ImGui_ImplGlfw_InitForOpenGL(pointer.window, true); ImGui::GetIO().IniFilename = nullptr;

    {
        // enter the render loop
        while (!glfwWindowShouldClose(pointer.window)) {

            // clear the color and depth buffer
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // render stuff
            gui();
            
            // swap buffers and poll events
            glfwSwapBuffers(pointer.window); glfwPollEvents();
        }
    }

    // destroy the GUI
    ImGui_ImplGlfw_Shutdown(); ImGui_ImplOpenGL3_Shutdown(); ImPlot::DestroyContext(); ImGui::DestroyContext();
}

void Viewer::gui() const {
    // begin the frame
    ImGui_ImplOpenGL3_NewFrame(); ImGui_ImplGlfw_NewFrame(); ImGui::NewFrame();

    // show demo
    ImGui::ShowDemoWindow(); ImPlot::ShowDemoWindow();

    // end the frame
    ImGui::Render(); ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
