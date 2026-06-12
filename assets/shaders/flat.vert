#version 410 core

layout(location = 0) in  vec3 a_Position;

layout(location = 0) out vec3 v_Position;

uniform mat4  u_Projection;
uniform mat4  u_View;
uniform mat4  u_Model;
uniform vec3  u_Color;

void main() {
    v_Position  = a_Position;
    gl_Position = u_Projection * u_View * u_Model * vec4(v_Position, 1.);
}
