#include "assignment5.hpp"
#include "node.hpp"
#include "parametric_shapes.hpp"

#include "config.hpp"
#include "external/glad/glad.h"
#include "core/Bonobo.h"
#include "core/FPSCamera.h"
#include "core/InputHandler.h"
#include "core/Log.h"
#include "core/LogView.h"
#include "core/Misc.h"
#include "core/utils.h"
#include "core/Window.h"
#include <imgui.h>
#include "external/imgui_impl_glfw_gl3.h"

#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "external/glad/glad.h"
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <cstdlib>
#include <stdexcept>
#include <stack>
#include <iostream>
#include <vector>
#include <stdlib.h>


enum class polygon_mode_t : unsigned int {
	fill = 0u,
	line,
	point
};

static polygon_mode_t get_next_mode(polygon_mode_t mode)
{
	return static_cast<polygon_mode_t>((static_cast<unsigned int>(mode) + 1u) % 3u);
}

eda221::Assignment5::Assignment5()
{
	Log::View::Init();

	window = Window::Create("EDA221: Assignment 5", config::resolution_x,
		config::resolution_y, config::msaa_rate, false);
	if (window == nullptr) {
		Log::View::Destroy();
		throw std::runtime_error("Failed to get a window: aborting!");
	}
	inputHandler = new InputHandler();
	window->SetInputHandler(inputHandler);
}

bool testCollision(glm::vec3 p1, float r1, glm::vec3 p2, float r2)
{
	if (glm::length(p1 - p2) < r1 + r2) {
		return true;
	}
	else {
		return false;
	}
}

bool testRaySphere(glm::vec3 pv, glm::vec3 v, glm::vec3 ps, float r)
{
	glm::vec3 rejection = (ps - pv) - ((pv * glm::dot(ps - pv, v)));
	if (glm::length(rejection) < r) {
		return true;
	}
	else return false;
}

bool visible(Node asteroid, Node camera,float tol) {
	if (asteroid.get_translation().z < camera.get_translation().z && asteroid.get_translation().z > camera.get_translation().z - tol) {
		return true;
	}
	else {
		return false;
	}
}

eda221::Assignment5::~Assignment5()
{
	delete inputHandler;
	inputHandler = nullptr;

	Window::Destroy(window);
	window = nullptr;

	Log::View::Destroy();
}

void
eda221::Assignment5::run()
{
	// Load the asteroid geometry
	auto asteroid_shape = parametric_shapes::createSphere(50u, 50u, 3.0f);
	if (asteroid_shape.vao == 0u) {
		LogError("Failed to retrieve the asteroid mesh");
		return;
	}

	// Load the world sphere/ the map geometry
	auto sphere_big_shape = parametric_shapes::createSphere(50u, 50u, 500.0f);
	if (sphere_big_shape.vao == 0u) {
		LogError("Failed to retrieve the big sphere mesh");
		return;
	}

	// Load bullet geometry
	auto bullet_shape = parametric_shapes::createSphere(10u, 10u, 2.0f);
	if (bullet_shape.vao == 0u) {
		LogError("Failed to retrieve the bullet mesh");
		return;
	}

	// Load the spaceship geometry
	auto const spaceship_shape = eda221::loadObjects("spaceship.obj");
	if (spaceship_shape.empty()) {
		LogError("Failed to load spaceship model");
		return;
	}

	// Set up the camera
	FPSCameraf mCamera(bonobo::pi / 4.0f,
		static_cast<float>(config::resolution_x) / static_cast<float>(config::resolution_y),
		0.01f, 1000.0f);
	mCamera.mWorld.SetTranslate(glm::vec3(0.0f, 0.0f, 15.0f));
	mCamera.mMouseSensitivity = 0.003f;
	mCamera.mMovementSpeed = 0.025;
	window->SetCamera(&mCamera);


	// Create the shader programs
	auto fallback_shader = eda221::createProgram("fallback.vert", "fallback.frag");
	if (fallback_shader == 0u) {
		LogError("Failed to load fallback shader");
		return;
	}

	auto cubemap_shader= eda221::createProgram("cubemap.vert", "cubemap.frag");
	if (cubemap_shader == 0u) {
		LogError("Failed to load cubemap shader");
		return;
	}

	auto spaceship_shader = eda221::createProgram("default.vert", "default.frag");
	if (spaceship_shader == 0u) {
		LogError("Failed to load default shader");
		return;
	}

	auto const light_position = glm::vec3(-2.0f, 4.0f, 2.0f);
	auto const set_uniforms = [&light_position](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
	};

	auto polygon_mode = polygon_mode_t::fill;

	// The spaceship
	auto spaceship = Node();
	spaceship.set_geometry(spaceship_shape.at(0));
	spaceship.set_program(spaceship_shader, set_uniforms);

	// A vector containing the asteroids
	unsigned int nbrofAsteroids = 10; // Number of asteroids
	auto asteroids = std::vector<Node>(nbrofAsteroids);
	// Loop filling the vector with asteroids
	for (int i = 0; i < asteroids.size(); i++) {
		asteroids[i] = Node();
		asteroids[i].set_geometry(asteroid_shape);
		asteroids[i].set_program(fallback_shader, set_uniforms);
	}

	// A bullet
	auto bullet = Node();
	bullet.set_geometry(bullet_shape);
	bullet.set_program(fallback_shader, set_uniforms);

	// The world sphere
	auto worldSphere = Node();
	worldSphere.set_geometry(sphere_big_shape);
	worldSphere.set_program(cubemap_shader, set_uniforms);

	auto sphere_cubemapTexture = loadTextureCubeMap("lightblue/posx.png", "lightblue/negx.png", "lightblue/posy.png", "lightblue/negy.png", "lightblue/negz.png", "lightblue/posz.png", true);
	auto spaceship_texture = loadTexture2D("metal-crate.png");
	worldSphere.add_texture("cubemap texture", sphere_cubemapTexture, GL_TEXTURE_CUBE_MAP);
	spaceship.add_texture("spaceship texture", spaceship_texture);

	// Camera node
	auto cameraNode = Node();
	cameraNode.add_child(&spaceship);

	glEnable(GL_DEPTH_TEST);

	// Enable face culling to improve performance
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_FRONT);
	//glCullFace(GL_BACK);

	f64 ddeltatime;
	size_t fpsSamples = 0;
	double nowTime, lastTime = GetTimeMilliseconds();
	double fpsNextTick = lastTime + 1000.0;

	
	glm::vec3 const cameraDirection = mCamera.mWorld.GetFront();
	bullet.translate(glm::vec3(0.0f,0.0f,-10.f));
	for (int i = 0; i < asteroids.size(); i++) {
			asteroids[i].set_translation(glm::vec3(rand() % 50 - 25.0f, rand() % 50 - 25.0f, -(rand() % 50 + 50.0f)));
	}

	while (!glfwWindowShouldClose(window->GetGLFW_Window())) {
		nowTime = GetTimeMilliseconds();
		ddeltatime = nowTime - lastTime;
		if (nowTime > fpsNextTick) {
			fpsNextTick += 1000.0;
			fpsSamples = 0;
		}
		fpsSamples++;

		glfwPollEvents();
		inputHandler->Advance();
		mCamera.Update(ddeltatime, *inputHandler);

		ImGui_ImplGlfwGL3_NewFrame();

		if (inputHandler->GetKeycodeState(GLFW_KEY_1) & JUST_PRESSED) {
			spaceship.set_program(fallback_shader, set_uniforms);
			worldSphere.set_program(cubemap_shader, set_uniforms);

		}
		if (inputHandler->GetKeycodeState(GLFW_KEY_Z) & JUST_PRESSED) {
			polygon_mode = get_next_mode(polygon_mode);
		}
		switch (polygon_mode) {
		case polygon_mode_t::fill:
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			break;
		case polygon_mode_t::line:
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			break;
		case polygon_mode_t::point:
			glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
			break;
		}

		cameraNode.set_translation(mCamera.mWorld.GetTranslation());
		spaceship.set_translation(glm::vec3(0.0f, 0.0f, -10.0f));
		bullet.set_translation(bullet.get_translation() + 0.01f * cameraDirection);
		//worldSphere.set_translation(cameraNode.get_translation() - glm::vec3(0.0f, 0.0f, nowTime/1000));

		for (int i = 0; i < asteroids.size(); i++) {
			if (!visible(asteroids[i], cameraNode, 100.0f)) {
				asteroids[i].set_translation(glm::vec3(cameraNode.get_translation().x + rand() % 50 - 25.0f, cameraNode.get_translation().y + rand() % 50 - 25.0f, cameraNode.get_translation().z - (rand() % 50 + 50.0f)));
			}
		}

		auto const window_size = window->GetDimensions();
		glViewport(0, 0, window_size.x, window_size.y);
		glClearDepthf(1.0f);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

		for (int i = 0; i < asteroids.size(); i++) {
			if (visible(asteroids[i], cameraNode, 200.0f)) {
				asteroids[i].render(mCamera.GetWorldToClipMatrix(), asteroids[i].get_transform());
			}
		}

		worldSphere.render(mCamera.GetWorldToClipMatrix(), worldSphere.get_transform());
		//bullet.render(mCamera.GetWorldToClipMatrix(), bullet.get_transform());

		// Traverse the scene graph and render all the nodes
		auto node_stack = std::stack<Node const*>();
		auto matrix_stack = std::stack<glm::mat4>();
		node_stack.push(&cameraNode);
		matrix_stack.push(glm::mat4());
		do {
			auto const* const current_node = node_stack.top();
			node_stack.pop();

			auto const parent_matrix = matrix_stack.top();
			matrix_stack.pop();

			auto const current_node_matrix = current_node->get_transform();

			//
			// Todo: Compute the current node's world matrix
			//
			auto const current_node_world_matrix = parent_matrix * current_node_matrix;
			current_node->render(mCamera.GetWorldToClipMatrix(), current_node_world_matrix);
			for (int i = static_cast<int>(current_node->get_children_nb()) - 1; i >= 0; --i) {
				node_stack.push(current_node->get_child(static_cast<size_t>(i)));
				matrix_stack.push(current_node_world_matrix);
			}
		} while (!node_stack.empty());

		/*if (testCollision(sphere.get_translation(), 3, cameraNode.get_translation() + spaceship.get_translation(), 2.1f)) {
			std::cout << "Crashed!";
		}*/

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		Log::View::Render();
		ImGui::Render();

		window->Swap();
		lastTime = nowTime;
	}

	glDeleteProgram(fallback_shader);
	fallback_shader = 0u;
}


int main()
{
	//! \todo Implement assignment 5
	Bonobo::Init();
	try {
		eda221::Assignment5 assignment5;
		assignment5.run();
	}
	catch (std::runtime_error const& e) {
		LogError(e.what());
	}
	Bonobo::Destroy();
}
