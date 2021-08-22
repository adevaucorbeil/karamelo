#pragma once

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <vector>

using namespace std;

void framebuffer_size_callback(GLFWwindow *window,
                               int width,
                               int height);

const char *vertexShaderSource = "#version 330 core\n"
                                 "layout (location = 0) in vec2 aPos;"
                                 "layout (location = 1) in vec3 aColor;"
                                 "out vec3 vertexColor;"
                                 "void main()"
                                 "{"
                                 "   gl_Position = vec4(aPos.x, aPos.y, 0.0, 1.0);"
                                     "vertexColor = aColor;"
                                 "}";
const char *fragmentShaderSource = "#version 330 core\n"
                                   "in vec3 vertexColor;"
                                   "out vec4 FragColor;"
                                   "void main()"
                                   "{"
                                   "   FragColor = vec4(vertexColor, 1.0f);"
                                   "}";

template<typename F>
int draw(int n, const float *const vertices, F update_vertices)
{
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  GLFWwindow *window = glfwCreateWindow(800, 600, "MPM", NULL, NULL);
  if (!window)
  {
    cout << "Failed to create GLFW window" << endl;
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    cout << "Failed to initialize GLAD" << endl;
    return -1;
  }

  unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
  glCompileShader(vertexShader);
  int success;
  char infoLog[512];
  glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
  if (!success)
  {
    glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
    cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << endl;
  }
  unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
  glCompileShader(fragmentShader);
  glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
  if (!success)
  {
    glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
    cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << endl;
  }
  unsigned int shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glLinkProgram(shaderProgram);
  glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
  if (!success) {
    glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
    cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << endl;
  }
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  unsigned int VAO, VBO;
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);
  glBindVertexArray(VAO);

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, 5*n*sizeof(float), NULL, GL_STREAM_DRAW);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void *)0);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void *)(2*sizeof(float)));
  glEnableVertexAttribArray(1);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  glPointSize(3.0f);

  while (!glfwWindowShouldClose(window))
  {
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    update_vertices(glfwGetKey(window, GLFW_KEY_R)                      == GLFW_PRESS,
                    glfwGetKey(window, GLFW_KEY_LEFT)                   == GLFW_PRESS,
                    glfwGetKey(window, GLFW_KEY_RIGHT)                  == GLFW_PRESS,
                    glfwGetKey(window, GLFW_KEY_UP)                     == GLFW_PRESS,
                    glfwGetKey(window, GLFW_KEY_DOWN)                   == GLFW_PRESS,
                    glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)  == GLFW_PRESS,
                    glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS,
                    xpos, ypos);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 5*n*sizeof(float), vertices);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glUseProgram(shaderProgram);
    glBindVertexArray(VAO);
    glDrawArrays(GL_POINTS, 0, n);

    glfwSwapBuffers(window);
    glfwPollEvents();
  }

  glfwTerminate();

  return 0;
}

void framebuffer_size_callback(GLFWwindow *window,
                               int width,
                               int height)
{
  glViewport(0, 0, width, height);
}