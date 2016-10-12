#include "assignment5.hpp"
#include "interpolation.hpp"
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

	window = Window::Create("EDA221: Assignment 3", config::resolution_x,
		config::resolution_y, config::msaa_rate, false);
	if (window == nullptr) {
		Log::View::Destroy();
		throw std::runtime_error("Failed to get a window: aborting!");
	}
	inputHandler = new InputHandler();
	window->SetInputHandler(inputHandler);
}

eda221::Assignment5::~Assignment5()
{
	delete inputHandler;
	inputHandler = nullptr;

	Window::Destroy(window);
	window = nullptr;
	Log::View::Destroy();
}

bool testSphereSphere(glm::vec3 p1, float r1, glm::vec3 p2, float r2)
{
	if (glm::length(p1 - p2) < (r1 + r2)) {
		return true;
	}
	return false;
}

bool testRaySphere(glm::vec3 pv, glm::vec3 ps, glm::vec3 v, float r1, float r2) {
	glm::vec3 rejection = ps - pv - v * (glm::dot(ps - pv, v));
	if (glm::length(rejection) < r1 + r2) {
		return true;
	}
	return false;
}

bool visible(Node asteroid, Node camera, float tol) {
	if (asteroid.get_translation().z < camera.get_translation().z && asteroid.get_translation().z > camera.get_translation().z - tol) {
		return true;
	}
	else {
		return false;
	}
}



void
eda221::Assignment5::run()
{
	// Load the sphere geometry
	auto circle_ring_shape = parametric_shapes::createCircleRing(4u, 60u, 1.0f, 2.0f);
	if (circle_ring_shape.vao == 0u) {
		LogError("Failed to retrieve the circle ring mesh");
		return;
	}

	float radius_ship = 0.2f;
	auto const sphere = parametric_shapes::createSphere(600u, 400u, radius_ship);
	if (sphere.vao == 0u)
		return;
	auto const sky_sphere = parametric_shapes::createSphere(600u, 400u, 300.0f);
	if (sphere.vao == 0u)
		return;

	//Load asteroid genometry
	float asteroid_radius = 1.0f;
	auto const asteroid_shape = parametric_shapes::createSphere(100u, 80u, asteroid_radius);
	if (sphere.vao == 0u)
		return;

	//Load bullet geometry
	float bullet_radius = 0.05f;
	auto const bullet_gem = parametric_shapes::createSphere(40u, 30u, bullet_radius);
	if (sphere.vao == 0u)
		return;

	// Set up the camera
	FPSCameraf mCamera(bonobo::pi / 4.0f,
		static_cast<float>(config::resolution_x) / static_cast<float>(config::resolution_y),
		0.01f, 1000.0f);
	mCamera.mWorld.SetTranslate(glm::vec3(0.0f, 0.0f, 100.0f));
	mCamera.mMouseSensitivity = 0.003f;
	mCamera.mMovementSpeed = 0.005;
	window->SetCamera(&mCamera);

	// Create the shader programs
	auto fallback_shader = eda221::createProgram("fallback.vert", "fallback.frag");
	if (fallback_shader == 0u) {
		LogError("Failed to load fallback shader");
		return;
	}

	auto shader = createProgram("default.vert", "default.frag");
	if (shader == 0u) {
		LogError("Failed to load shader");
		return;
	}

	GLuint diffuse_shader = 0u, normal_shader = 0u, texcoord_shader = 0u, phong_shader = 0u, bumpmap_shader = 0u, cubemap_shader = 0u;
	auto const reload_shaders = [&diffuse_shader, &normal_shader, &texcoord_shader, &phong_shader, &bumpmap_shader, &cubemap_shader]() {
		if (diffuse_shader != 0u)
			glDeleteProgram(diffuse_shader);
		diffuse_shader = eda221::createProgram("diffuse.vert", "diffuse.frag");
		if (diffuse_shader == 0u)
			LogError("Failed to load diffuse shader");

		if (normal_shader != 0u)
			glDeleteProgram(normal_shader);
		normal_shader = eda221::createProgram("normal.vert", "normal.frag");
		if (normal_shader == 0u)
			LogError("Failed to load normal shader");

		if (texcoord_shader != 0u)
			glDeleteProgram(texcoord_shader);
		texcoord_shader = eda221::createProgram("texcoord.vert", "texcoord.frag");
		if (texcoord_shader == 0u)
			LogError("Failed to load texcoord shader");
		phong_shader = eda221::createProgram("phong.vert", "phong.frag");
		if (phong_shader == 0u) {
			LogError("Failed to load phong shader");
		}
		bumpmap_shader = eda221::createProgram("bumpmap.vert", "bumpmap.frag");
		if (bumpmap_shader == 0u) {
			LogError("Failed to load bumpmap shader");
		}
		cubemap_shader = eda221::createProgram("cubemap.vert", "cubemap.frag");
		if (cubemap_shader == 0u) {
			LogError("Failed to load cubemap shader");
		}
	};
	reload_shaders();

	auto light_position = glm::vec3(-20.0f, 40.0f, 20.0f);
	auto const set_uniforms = [&light_position](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
	};

	auto camera_position = mCamera.mWorld.GetTranslation();
	auto ambient = glm::vec3(0.1f, 0.1f, 0.1f);
	auto diffuse = glm::vec3(0.0f, 0.2f, 0.5f);
	auto specular = glm::vec3(1.0f, 1.0f, 1.0f);
	auto shininess = 50.0f;
	auto const phong_set_uniforms = [&light_position, &camera_position, &ambient, &diffuse, &specular, &shininess](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform3fv(glGetUniformLocation(program, "CamPos"), 1, glm::value_ptr(camera_position));
		glUniform3fv(glGetUniformLocation(program, "ka"), 1, glm::value_ptr(ambient));
		glUniform3fv(glGetUniformLocation(program, "kd"), 1, glm::value_ptr(diffuse));
		glUniform3fv(glGetUniformLocation(program, "ks"), 1, glm::value_ptr(specular));
		glUniform1f(glGetUniformLocation(program, "shininess"), shininess);
	};

	auto diffuse_asteroid = glm::vec3(0.3f, 0.2f, 0.2f);
	auto const asteroid_uniforms = [&light_position, &camera_position, &ambient, &diffuse_asteroid, &specular, &shininess](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform3fv(glGetUniformLocation(program, "CamPos"), 1, glm::value_ptr(camera_position));
		glUniform3fv(glGetUniformLocation(program, "ka"), 1, glm::value_ptr(ambient));
		glUniform3fv(glGetUniformLocation(program, "kd"), 1, glm::value_ptr(diffuse_asteroid));
		glUniform3fv(glGetUniformLocation(program, "ks"), 1, glm::value_ptr(specular));
		glUniform1f(glGetUniformLocation(program, "shininess"), shininess);
	};

	auto diffuse_bullet = glm::vec3(1.0f, 0.0f, 0.0f);
	auto const bullet_uniforms = [&light_position, &camera_position, &ambient, &diffuse_bullet, &specular, &shininess](GLuint program) {
		glUniform3fv(glGetUniformLocation(program, "light_position"), 1, glm::value_ptr(light_position));
		glUniform3fv(glGetUniformLocation(program, "CamPos"), 1, glm::value_ptr(camera_position));
		glUniform3fv(glGetUniformLocation(program, "ka"), 1, glm::value_ptr(ambient));
		glUniform3fv(glGetUniformLocation(program, "kd"), 1, glm::value_ptr(diffuse_bullet));
		glUniform3fv(glGetUniformLocation(program, "ks"), 1, glm::value_ptr(specular));
		glUniform1f(glGetUniformLocation(program, "shininess"), shininess);
	};

	int spaceship_has_tex = 1;
	auto const spaceship_uniforms = [&spaceship_has_tex](GLuint program) {
		glUniform1i(glGetUniformLocation(program, "has_texture"), spaceship_has_tex);
	};

	auto polygon_mode = polygon_mode_t::fill;

	// A vector containing the asteroids
	unsigned int nbrofAsteroids = 40; // Number of asteroids
	auto asteroids = std::vector<Node>(nbrofAsteroids);
	auto asteroids_alive = std::vector<bool>(nbrofAsteroids);
	auto asteroid_bump = loadTexture2D("stone47_bump.png");
	// Loop filling the vector with asteroids
	for (int i = 0; i < asteroids.size(); i++) {
		asteroids[i] = Node();
		asteroids[i].set_geometry(asteroid_shape);
		//asteroids[i].set_program(fallback_shader, set_uniforms);
		asteroids[i].set_program(bumpmap_shader, asteroid_uniforms);
		asteroids[i].add_texture("asteroid_bump", asteroid_bump);
		asteroids_alive[i] = false;
	}

	//Create spaceship
	auto const spaceship_shape = eda221::loadObjects("spaceship.obj");
	if (spaceship_shape.empty()) {
		LogError("Failed to load spaceship model");
		return;
	}
	auto spaceship = Node();
	spaceship.set_geometry(spaceship_shape.at(0));
	spaceship.set_program(shader, spaceship_uniforms);
	auto spaceship_tex = loadTexture2D("metal-crate.png");
	spaceship.add_texture("metal", spaceship_tex);
	//spaceship.rotate_x(1.0f);	

	auto sphere_node = Node();
	sphere_node.set_geometry(sphere);
	sphere_node.set_program(fallback_shader, set_uniforms);

	auto sphere_tex = loadTexture2D("earth_bump.png");
	sphere_node.add_texture("sphere tex", sphere_tex);

	//Creates a camera node to attach spaceship to
	auto camera_node = Node();
	camera_node.add_child(&spaceship);
	//camera_node.set_translation(mCamera.mWorld.GetTranslation());

	//skybox
	auto skybox_node = Node();
	skybox_node.set_geometry(sky_sphere);
	skybox_node.set_program(cubemap_shader, set_uniforms);

	auto skybox = loadTextureCubeMap("snow/posx.png", "snow/negx.png", "snow/posy.png", "snow/negy.png", "snow/negz.png", "snow/posz.png", true);
	//auto skybox = loadTextureCubeMap("sunset_sky/posx.png", "sunset_sky/negx.png", "sunset_skynow/posy.png", "sunset_sky/negy.png", "sunset_sky/negz.png", "sunset_sky/posz.png", true);
	skybox_node.add_texture("skybox", skybox, GL_TEXTURE_CUBE_MAP);

	//Set up bullets
	int nbrofBullets = 4;
	auto bullets = std::vector<Node>(nbrofBullets);
	auto bullets_alive = std::vector<bool>(nbrofBullets);
	auto bullet_hit = std::vector<bool>(nbrofBullets);
	for (int k = 0; k < nbrofBullets; k++) {
		bullets[k].set_geometry(bullet_gem);
		bullets[k].set_program(fallback_shader, set_uniforms);
		bullets_alive[k] = false;
		bullet_hit[k] = false;
	}

	//glm::vec3 direction = glm::vec3(0.0f);
	//glm::vec3 startpoint = glm::vec3(0.0f);
	////bool render_bullet = false;


	//Set up bullet
	auto bullet = Node();
	bullet.set_geometry(bullet_gem);
	bullet.set_program(fallback_shader, set_uniforms);
	glm::vec3 direction = glm::vec3(0.0f);
	glm::vec3 startpoint = glm::vec3(0.0f);
	//bool render_bullet = false;
	bool bullet_alive = false;
	bool hit;


	glEnable(GL_DEPTH_TEST);
	// Enable face culling to improve performance:
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_FRONT);
	//glCullFace(GL_BACK);


	f64 ddeltatime;
	size_t fpsSamples = 0;
	double nowTime, lastTime = GetTimeMilliseconds();
	double fpsNextTick = lastTime + 1000.0;

	int i = 0;
	glm::vec3 const cam = mCamera.mWorld.GetFront();
	bullet.translate(glm::vec3(0.0, 0.0, 4.0));
	auto controlPoints = std::vector<glm::vec3>(nbrofBullets * 2);
	auto path_pos = std::vector<float>(nbrofBullets);
	for (int k = 0; k < nbrofBullets; k++) {
		path_pos[k] = 0.0f;
	}
	float path_velocity = 0.01f;
	int bullet_counter = 0;
	auto interpol = std::vector<glm::vec3>(nbrofBullets);

	bool collision = false;

	glm::vec3 spaceship_offset = glm::vec3(0.0, -3.0, -10.0);
	//added speed
	float dt = 0.0f;
	int score = 0;

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

		//sphere_node.set_program(bumpmap_shader, phong_set_uniforms);
		
		//sphere_node.set_translation(glm::vec3(0.0, -2.0, -5.0));
		//asteroid.set_translation(glm::vec3(0.0, 0.0, -3.0));
		//bullet.set_program(phong_shader, bullet_uniforms);

		spaceship.set_translation(spaceship_offset);
		spaceship.set_scaling(glm::vec3(0.5f));

		/*bool collision = testSphereSphere(sphere_node.get_translation(), sphere_node.get_scaling()[1] * radius_ship,
		asteroid.get_translation(), asteroid.get_scaling[1] * asteroid_radius);*/
		//Check collision with all asteroids
		for (int k = 0; k < nbrofAsteroids; k++) {
			collision = testSphereSphere(camera_node.get_translation() + spaceship.get_translation(), 2.1f * spaceship.get_scaling().x,
				asteroids[k].get_translation(), asteroid_radius * asteroids[k].get_scaling().x);
			if (collision && asteroids_alive[k]) {
				asteroids_alive[k] = false;
				printf("COLLISION!!!  ");
				score -= 1;
			}
		}

		//Shoot bullets
		bool shoot_fired = false;
		for (int k = 0; k < nbrofBullets; k++) {
			if (i % 10 == 0 && !bullets_alive[k] && !shoot_fired) {
				controlPoints[bullet_counter] = camera_node.get_translation() + spaceship.get_translation();
				controlPoints[bullet_counter + 1] = controlPoints[bullet_counter] - (spaceship_offset - glm::vec3(0.0f, 0.0f, spaceship_offset.z)) 
					+ mCamera.mWorld.GetFront() * 20.0f;
				bullets_alive[k] = true;
				shoot_fired = true;
				bullet_counter += 2;
			}
			//if (i % 100 == 0) {
			//	controlPoints[0] = mCamera.mWorld.GetTranslation() + spaceship.get_translation();
			//	controlPoints[1] = controlPoints[0] - (spaceship_offset - glm::vec3(0.0f, 0.0f, spaceship_offset.z))
			//		+ mCamera.mWorld.GetFront() * 20.0f;
			//	bullets_alive[k] = true;
			//	//bullet_counter += 2;
			//}
			//if (i % 100 == 25) {
			//	controlPoints[2] = mCamera.mWorld.GetTranslation() + spaceship.get_translation();
			//	controlPoints[3] = controlPoints[2] - (spaceship_offset - glm::vec3(0.0f, 0.0f, spaceship_offset.z))
			//		+ mCamera.mWorld.GetFront() * 20.0f;
			//	bullets_alive[k] = true;
			//	//bullet_counter += 2;
			//}
			if (bullet_counter >= nbrofBullets * 2) {
				bullet_counter = 0;
			}
		}
		shoot_fired = false;

		//Shoot bullet
		//if ((int)nowTime % 1000 < 0.01) {
		//if (i % 100 == 0) {
		//	controlPoints[0] = mCamera.mWorld.GetTranslation() + spaceship.get_translation();
		//	controlPoints[1] = controlPoints[0] - (spaceship_offset - glm::vec3(0.0f, 0.0f, spaceship_offset.z)) + mCamera.mWorld.GetFront() * 20.0f;
		//	bullet_alive = true;
		//}

		//glm::vec3 interpol;
		for (int k = 0; k < nbrofBullets; k++) {
			if (bullets_alive[k]) {
				int j = floor(path_pos[k]);
				if (j >= 1) {
					bullets_alive[k] = false;
				}
				interpol[k] = interpolation::evalLERP(controlPoints[k * 2], controlPoints[k * 2 + 1], path_pos[k] - j);
				bullets[k].set_translation(interpol[k]);
				path_pos[k] += path_velocity;
			}
		}



		//Check bullet hit
		for (int k = 0; k < nbrofAsteroids; k++) {
			for (int j = 0; j < nbrofBullets; j++) {
				bullet_hit[j] = testRaySphere(bullets[j].get_translation(), asteroids[k].get_translation(), 
					direction, asteroid_radius * asteroids[k].get_scaling().x, bullet_radius);
				if (bullet_hit[j] && bullets_alive[j] && asteroids_alive[k]) {
					bullets_alive[j] = false;
					printf("BOOM!  ");
					asteroids_alive[k] = false;
					score += 1;
				}
			}
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

		//camera_position = mCamera.mWorld.GetTranslation();
		mCamera.mWorld.Translate(glm::vec3(0.0f, 0.0f, dt));
		camera_node.set_translation(mCamera.mWorld.GetTranslation());
		//camera_node.set_rotation(mCamera.mWorld.GetFront());
		if (i % 100 == 0) {
			dt += -0.01f;
		}


		for (int k = 0; k < asteroids.size(); k++) {
			if (!visible(asteroids[k], camera_node, 100.0f) || !asteroids_alive[k]) {
				asteroids[k].set_translation(glm::vec3(camera_node.get_translation().x + rand() % 50 - 25.0f, camera_node.get_translation().y + rand() % 50 - 25.0f, camera_node.get_translation().z - (rand() % 50 + 50.0f)));
				asteroids_alive[k] = true;
			}
		}

		auto const window_size = window->GetDimensions();
		glViewport(0, 0, window_size.x, window_size.y);
		glClearDepthf(1.0f);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

		//sphere_node.render(mCamera.GetWorldToClipMatrix(), sphere_node.get_transform());
		skybox_node.render(mCamera.GetWorldToClipMatrix(), skybox_node.get_transform());
		for (int k = 0; k < nbrofAsteroids; k++) {
			if (visible(asteroids[k], camera_node, 200.0f) && asteroids_alive[k]) {
				asteroids[k].render(mCamera.GetWorldToClipMatrix(), asteroids[k].get_transform());
			}
		}
		for (int k = 0; k < nbrofBullets; k++){
			if (bullets_alive[k]) {
				bullets[k].render(mCamera.GetWorldToClipMatrix(), bullets[k].get_transform());
			}
		}
		

		// Traverse the scene graph and render all the nodes
		auto node_stack = std::stack<Node const*>();
		auto matrix_stack = std::stack<glm::mat4>();
		node_stack.push(&camera_node);
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


		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		Log::View::Render();

		bool opened = ImGui::Begin("Scene Control", &opened, ImVec2(300, 100), -1.0f, 0);
		if (opened) {
			ImGui::ColorEdit3("Ambient", glm::value_ptr(ambient));
			ImGui::ColorEdit3("Diffuse", glm::value_ptr(diffuse));
			ImGui::ColorEdit3("Specular", glm::value_ptr(specular));
			ImGui::SliderFloat("Shininess", &shininess, 0.0f, 1000.0f);
			ImGui::SliderFloat3("Light Position", glm::value_ptr(light_position), -20.0f, 20.0f);
			//			ImGui::SliderInt("Faces Nb", &faces_nb, 1u, 16u);
		}
		ImGui::End();

		ImGui::Begin("Render Time", &opened, ImVec2(120, 50), -1.0f, 0);
		if (opened)
			ImGui::Text("%.3f ms", ddeltatime);
		ImGui::End();
		ImGui::Begin("Score ", &opened, ImVec2(120, 50), -10.0f, 0);
		if (opened)
			ImGui::Text("Your score: %d", score);
		ImGui::End();

		ImGui::Render();

		window->Swap();
		lastTime = nowTime;

		i++;
	}
	printf("Your score: %d ", score);

	glDeleteProgram(texcoord_shader);
	texcoord_shader = 0u;
	glDeleteProgram(normal_shader);
	normal_shader = 0u;
	glDeleteProgram(diffuse_shader);
	diffuse_shader = 0u;
	glDeleteProgram(fallback_shader);
	diffuse_shader = 0u;
	glDeleteProgram(phong_shader);
	phong_shader = 0u;
	glDeleteProgram(bumpmap_shader);
	bumpmap_shader = 0u;
	glDeleteProgram(cubemap_shader);
	cubemap_shader = 0u;
}

int main()
{
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
