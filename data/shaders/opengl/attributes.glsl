// Copyright © 2008-2015 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

uniform mat4 uProjectionMatrix;
uniform mat4 uViewMatrix;
uniform mat4 uViewMatrixInverse;
uniform mat4 uViewProjectionMatrix;
uniform mat3 uNormalMatrix;

//Light uniform parameters
struct Light {
	vec4 diffuse;
	vec4 specular;
	vec4 position;
};
uniform Light uLight[4];

struct Material {
	vec4 emission;
	vec4 ambient;
	vec4 diffuse;
	vec4 specular;
	float shininess;
};

#ifdef VERTEX_SHADER

in vec4 a_vertex;
in vec3 a_normal;
in vec4 a_color;
in vec4 a_uv0;
in mat4 a_transform;

#endif // VERTEX_SHADER