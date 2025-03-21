#version 330 core
// please send comments or report bug to: liangliang.nan@gmail.com

// It uses geometry shader for vertex generation.
// The code does not cover round caps.

//#extension GL_EXT_gpu_shader4 : enable

uniform mat4    MV;

in  vec3 vtx_position;	// point position
in  vec2 vtx_texcoord;

out vec2 vOutTexcoord;

void main()
{
    gl_Position = MV * vec4(vtx_position, 1.0);

    vOutTexcoord = vtx_texcoord;

//	if (clippingPlaneEnabled) {
//		gl_ClipDistance[0] = dot(new_position, clippingPlane0);
//		if (crossSectionEnabled)
//			gl_ClipDistance[1] = dot(new_position, clippingPlane1);
//	}
}
