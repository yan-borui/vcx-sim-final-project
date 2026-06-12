#version 410 core

layout(location = 0) in  vec3 a_Position;
layout(location = 1) in  vec3 a_Normal;
layout(location = 2) in  vec3 a_Offset;
layout(location = 3) in  vec3 a_Color;

layout(location = 0) out vec3 v_Position;
layout(location = 1) out vec3 v_Normal;
layout(location = 2) out vec3 v_Color;

uniform mat4 u_FluidProjection;
uniform mat4 u_FluidView;

void main() {
    v_Position  = a_Position + a_Offset;
    v_Normal    = a_Normal;
    v_Color     = a_Color;
    gl_Position = u_FluidProjection * u_FluidView * vec4(v_Position, 1.);
}
