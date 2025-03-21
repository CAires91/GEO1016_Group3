/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <easy3d/viewer/viewer.h>

#include <thread>
#include <chrono>
#include <iostream>


#if !defined(_WIN32)
#  include <unistd.h>
#  include <sys/wait.h>
#endif


#include <easy3d/viewer/opengl.h>		// Initialize with glewInit()
#include <3rd_party/glfw/include/GLFW/glfw3.h>	// Include glfw3.h after our OpenGL definitions

#include <easy3d/core/surface_mesh.h>
#include <easy3d/core/graph.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/viewer/drawable_points.h>
#include <easy3d/viewer/drawable_lines.h>
#include <easy3d/viewer/drawable_triangles.h>
#include <easy3d/viewer/renderer.h>
#include <easy3d/viewer/shader_program.h>
#include <easy3d/viewer/shader_manager.h>
#include <easy3d/viewer/transform.h>
#include <easy3d/viewer/primitives.h>
#include <easy3d/viewer/camera.h>
#include <easy3d/viewer/manipulated_camera_frame.h>
#include <easy3d/viewer/key_frame_interpolator.h>
#include <easy3d/viewer/framebuffer_object.h>
#include <easy3d/viewer/opengl_error.h>
#include <easy3d/viewer/opengl_timer.h>
#include <easy3d/viewer/setting.h>
#include <easy3d/fileio/resources.h>
#include <easy3d/fileio/surface_mesh_io.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <easy3d/util/dialogs.h>
#include <easy3d/util/file_system.h>
#include <easy3d/util/logging.h>
#include <easy3d/util/timer.h>


 // enforce the same behavior on macOS and other platforms (i.e., Windows, Linux)
#ifdef __APPLE__
#define EASY3D_MOD_CONTROL GLFW_MOD_SUPER
#else
#define EASY3D_MOD_CONTROL GLFW_MOD_CONTROL
#endif


namespace easy3d {


	Viewer::Viewer(
		const std::string& title /* = "easy3d Viewer" */,
		int samples /* = 4 */,
		int gl_major /* = 3 */,
		int gl_minor /* = 2 */,
		bool full_screen /* = false */,
		bool resizable /* = true */,
		int depth_bits /* = 24 */,
		int stencil_bits /* = 8 */
	)
		: title_(title)
		, camera_(nullptr)
		, samples_(0)
		, full_screen_(full_screen)
		, process_events_(true)
		, gpu_timer_(nullptr)
		, button_(-1)
		, modifiers_(-1)
		, drag_active_(false)
		, mouse_x_(0)
		, mouse_y_(0)
		, mouse_pressed_x_(0)
		, mouse_pressed_y_(0)
		, show_pivot_point_(false)
		, drawable_axes_(nullptr)
		, show_camera_path_(false)
		, model_idx_(-1)
	{
		// Avoid locale-related number parsing issues.
		setlocale(LC_NUMERIC, "C");

		// create and setup window
		window_ = create_window(title, samples, gl_major, gl_minor, full_screen, resizable,
			depth_bits, stencil_bits);
		setup_callbacks(window_);

		// create and setup the camera
		camera_ = new Camera;
		camera_->setType(Camera::PERSPECTIVE);
		camera_->setUpVector(vec3(0, 1, 0)); // y pointing up
		camera_->setViewDirection(vec3(0, 0, -1)); // Z pointing out
		camera_->showEntireScene();
		camera_->connect(this, &Viewer::update);

		int fw, fh;
		glfwGetFramebufferSize(window_, &fw, &fh);
		// needs to be executed once to ensure the viewer is initialized with correct viewer size
		callback_event_resize(fw, fh);

		// create a GPU timer
		gpu_timer_ = new OpenGLTimer(false);

		// initialize various variables
		background_color_ = vec4(0.9f, 0.9f, 1.0f, 1.0f);

		mouse_x_ = mouse_y_ = 0;
		mouse_pressed_x_ = mouse_pressed_y_ = 0;
		button_ = -1;
		modifiers_ = 0;
		drag_active_ = false;
		process_events_ = true;

		/* Poll for events once before starting a potentially lengthy loading process.*/
		glfwPollEvents();
	}


	GLFWwindow* Viewer::create_window(
		const std::string& title,
		int samples,
		int gl_major,
		int gl_minor,
		bool full_screen,
		bool resizable,
		int depth_bits,
		int stencil_bits)
	{
		glfwSetErrorCallback(
			[](int error, const char *desc) {
			if (error == GLFW_NOT_INITIALIZED)
				LOG(ERROR) << "GLFW error " << error << ": " << desc;
		});

		if (!glfwInit()) {
			LOG(ERROR) << "could not initialize GLFW!";
			throw std::runtime_error("could not initialize GLFW!");
		}

		glfwSetTime(0);

		// Reset the hints, allowing viewers to have different hints.
		glfwDefaultWindowHints();

		glfwWindowHint(GLFW_SAMPLES, samples);

		glfwWindowHint(GLFW_STENCIL_BITS, stencil_bits);
		glfwWindowHint(GLFW_DEPTH_BITS, depth_bits);

		/* Request a forward compatible OpenGL glMajor.glMinor core profile context.
		   Default value is an OpenGL 3.2 core profile context. */
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, gl_major);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, gl_minor);

#ifdef __APPLE__
		// The only OpenGL 3.x and 4.x contexts currently supported by macOS are
		// forward-compatible, core profile contexts.
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#else
		if (gl_major >= 3) {
			if (gl_minor >= 2)
				glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);	// 3.2+ only
			if (gl_minor >= 0)
				glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);			// 3.0+ only
		}
#endif

		// make the whole window transparent
		//glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GLFW_TRUE);

		glfwWindowHint(GLFW_VISIBLE, GL_FALSE);
		glfwWindowHint(GLFW_RESIZABLE, resizable ? GL_TRUE : GL_FALSE);

		int w = 960;
		int h = 800;
		GLFWwindow*	window = nullptr;
		if (full_screen) {
			GLFWmonitor *monitor = glfwGetPrimaryMonitor();
			const GLFWvidmode *mode = glfwGetVideoMode(monitor);
			w = mode->width;
			h = mode->height;
			window = glfwCreateWindow(w, h, title.c_str(), monitor, nullptr);
		}
		else {
			window = glfwCreateWindow(w, h, title.c_str(), nullptr, nullptr);
		}

		if (!window) {
			glfwTerminate();
			LOG(ERROR) << "could not create an OpenGL " << std::to_string(gl_major)
				<< "." << std::to_string(gl_minor) << " context!";
			throw std::runtime_error("could not create an OpenGL " +
				std::to_string(gl_major) + "." + std::to_string(gl_minor) + " context!");
		}
		glfwSetWindowUserPointer(window, this);
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

		// Enable vsync
		glfwMakeContextCurrent(window);
		glfwSwapInterval(1);

		// Load OpenGL and its extensions
		if (glewInit() != GLEW_OK) {
			glGetError(); // pull and ignore unhandled errors like GL_INVALID_ENUM
			throw std::runtime_error("failed to load OpenGL and its extensions!");
		}

#ifndef NDEBUG
		opengl::setup_gl_debug_callback();
#endif

		LOG(INFO) << "OpenGL vendor:            " << glGetString(GL_VENDOR);
		LOG(INFO) << "OpenGL renderer:          " << glGetString(GL_RENDERER);
		LOG(INFO) << "OpenGL version requested: " << gl_major << "." << gl_minor;
		LOG(INFO) << "OpenGL version received:  " << glGetString(GL_VERSION);
		LOG(INFO) << "GLSL version received:    " << glGetString(GL_SHADING_LANGUAGE_VERSION);

		glGetIntegerv(GL_SAMPLES, &samples_);
		int max_num = 0;
		glGetIntegerv(GL_MAX_SAMPLES, &max_num);

		// warn the user if the requests were not satisfied
		if (samples > 0 && samples_ != samples) {
			if (samples_ == 0)
				LOG(WARNING) << "MSAA is not available (" << samples << " samples requested)";
			else
				LOG(WARNING) << "MSAA is available with " << samples_ << " samples (" << samples << " requested but max support is " << max_num << ")";
		}
		else
			LOG(INFO) << "samples received:         " << samples_ << " (" << samples << " requested, max support is " << max_num << ")";

		float xscale, yscale;
		glfwGetWindowContentScale(window, &xscale, &yscale);
		dpi_scaling_ = static_cast<double>(xscale + yscale) * 0.5;
		LOG(INFO) << "DPI scaling:              " << dpi_scaling();

		return window;
	}



	void Viewer::setup_callbacks(GLFWwindow *window) {
		glfwSetCursorPosCallback(window, [](GLFWwindow *win, double x, double y)
		{
			auto viewer = reinterpret_cast<Viewer*>(glfwGetWindowUserPointer(win));
			if (!viewer->process_events_)
				return;

			int w, h;
			glfwGetWindowSize(win, &w, &h);
			if (x >= 0 && x <= w && y >= 0 && y <= h)
				viewer->callback_event_cursor_pos(x, y);
			else if (viewer->drag_active_) {
				// Restrict the cursor to be within the client area during dragging
				if (x < 0) x = 0;	if (x > w) x = w;
				if (y < 0) y = 0;	if (y > h) y = h;
				glfwSetCursorPos(win, x, y);
			}
		});

		glfwSetMouseButtonCallback(window, [](GLFWwindow *win, int button, int action, int modifiers)
		{
			auto viewer = reinterpret_cast<Viewer*>(glfwGetWindowUserPointer(win));
			if (!viewer->process_events_)
				return;
			viewer->callback_event_mouse_button(button, action, modifiers);
		});

		glfwSetKeyCallback(window, [](GLFWwindow *win, int key, int scancode, int action, int mods)
		{
			auto viewer = reinterpret_cast<Viewer*>(glfwGetWindowUserPointer(win));
			if (!viewer->process_events_)
				return;
			(void)scancode;
			viewer->callback_event_keyboard(key, action, mods);
		});

		glfwSetCharCallback(window, [](GLFWwindow *win, unsigned int codepoint)
		{
			auto viewer = reinterpret_cast<Viewer*>(glfwGetWindowUserPointer(win));
			if (!viewer->process_events_)
				return;
			viewer->callback_event_character(codepoint);
		});

		glfwSetDropCallback(window, [](GLFWwindow *win, int count, const char **filenames)
		{
			auto viewer = reinterpret_cast<Viewer*>(glfwGetWindowUserPointer(win));
			if (!viewer->process_events_)
				return;
			viewer->callback_event_drop(count, filenames);
		});

		glfwSetScrollCallback(window, [](GLFWwindow *win, double dx, double dy)
		{
			auto viewer = reinterpret_cast<Viewer*>(glfwGetWindowUserPointer(win));
			if (!viewer->process_events_)
				return;
			viewer->callback_event_scroll(dx, dy);
		});

		/* React to framebuffer size events -- includes window size events and also
		 * catches things like dragging a window from a Retina-capable screen to a
		 * normal screen on Mac OS X */
		glfwSetFramebufferSizeCallback(window, [](GLFWwindow *win, int width, int height)
		{
			auto viewer = reinterpret_cast<Viewer*>(glfwGetWindowUserPointer(win));
			if (!viewer->process_events_)
				return;
			viewer->callback_event_resize(width, height);
		});

		// notify when the screen has lost focus (e.g. application switch)
		glfwSetWindowFocusCallback(window, [](GLFWwindow *win, int focused)
		{
			auto viewer = reinterpret_cast<Viewer*>(glfwGetWindowUserPointer(win));
			viewer->focus_event(focused != 0);// true for focused
		});

		glfwSetWindowCloseCallback(window, [](GLFWwindow *win)
		{
			glfwSetWindowShouldClose(win, true);
		});
	}


	bool Viewer::callback_event_cursor_pos(double x, double y) {
		int px = static_cast<int>(x);
		int py = static_cast<int>(y);
		try {
			int dx = px - mouse_x_;
			int dy = py - mouse_y_;
			mouse_x_ = px;
			mouse_y_ = py;
			if (drag_active_)
				return mouse_drag_event(px, py, dx, dy, button_, modifiers_);
			else
				return mouse_free_move_event(px, py, dx, dy, modifiers_);
		}
		catch (const std::exception &e) {
			LOG(ERROR) << "Caught exception in event handler: " << e.what();
			return false;
		}
	}


	bool Viewer::callback_event_mouse_button(int button, int action, int modifiers) {
		try {
			if (action == GLFW_PRESS) {
				drag_active_ = true;
				button_ = button;
				modifiers_ = modifiers;
				mouse_pressed_x_ = mouse_x_;
				mouse_pressed_y_ = mouse_y_;
				return mouse_press_event(mouse_x_, mouse_y_, button, modifiers);
			}
			else if (action == GLFW_RELEASE) {
				drag_active_ = false;
				return mouse_release_event(mouse_x_, mouse_y_, button, modifiers);
			}
			else {
				drag_active_ = false;
				LOG(INFO) << "GLFW_REPEAT? Seems never happen";
				return false;
			}
		}
		catch (const std::exception &e) {
			LOG(ERROR) << "Caught exception in event handler: " << e.what();
			return false;
		}
	}


	bool Viewer::callback_event_keyboard(int key, int action, int modifiers) {
		try {
			if (action == GLFW_PRESS || action == GLFW_REPEAT) {
				return key_press_event(key, modifiers);
			}
			else {
				return key_release_event(key, modifiers);
			}
		}
		catch (const std::exception &e) {
			LOG(ERROR) << "Caught exception in event handler: " << e.what();
			return false;
		}
	}


	bool Viewer::callback_event_character(unsigned int codepoint) {
		try {
			return char_input_event(codepoint);
		}
		catch (const std::exception &e) {
			LOG(ERROR) << "Caught exception in event handler: " << e.what();
			return false;
		}
	}


	bool Viewer::callback_event_drop(int count, const char **filenames) {
		try {
			std::vector<std::string> arg(count);
			for (int i = 0; i < count; ++i)
				arg[i] = filenames[i];
			return drop_event(arg);
		}
		catch (const std::exception &e) {
			LOG(ERROR) << "Caught exception in event handler: " << e.what();
			return false;
		}
	}


	bool Viewer::callback_event_scroll(double dx, double dy) {
		try {
			return mouse_scroll_event(mouse_x_, mouse_y_, static_cast<int>(dx), static_cast<int>(dy));
		}
		catch (const std::exception &e) {
			LOG(ERROR) << "Caught exception in event handler: " << e.what();
			return false;
		}
	}


	void Viewer::callback_event_resize(int w, int h) {
		if (w == 0 && h == 0)
			return;

		try {
			// camera is manipulated by the mouse, working in the screen coordinate system
			int win_w, win_h;
			glfwGetWindowSize(window_, &win_w, &win_h);
			camera_->setScreenWidthAndHeight(win_w, win_h);
			glViewport(0, 0, w, h);
			post_resize(w, h);
		}
		catch (const std::exception &e) {
			LOG(ERROR) << "Caught exception in event handler: " << e.what();
		}
	}


	bool Viewer::focus_event(bool focused) {
		if (focused) {
			// ...
		}
		return false;
	}


	Viewer::~Viewer() {
		cleanup();
        LOG(INFO) << "viewer terminated. Bye!";
	}

	void Viewer::clear_scene() {
		for (auto m : models_)
			delete m;
		models_.clear();

		for (auto d : drawables_)
			delete d;
		drawables_.clear();
	}

	void Viewer::cleanup() {
		// viewer may have already been destroyed by the user
		if (!window_)
			return;

		if (camera_) { delete camera_; camera_ = nullptr; }
		if (drawable_axes_) { delete drawable_axes_; drawable_axes_ = nullptr; }

		clear_scene();

		if (gpu_timer_) { delete gpu_timer_; gpu_timer_ = nullptr; }

		ShaderManager::terminate();

		glfwDestroyWindow(window_);
		window_ = nullptr;
		glfwTerminate();
	}


	void Viewer::set_title(const std::string &title) {
		if (title != title_) {
			glfwSetWindowTitle(window_, title.c_str());
			title_ = title;
		}
	}


	void Viewer::resize(int w, int h) {
#ifndef __APPLE__
		w = static_cast<int>(w * dpi_scaling());
		h = static_cast<int>(h * dpi_scaling());
#endif
		glfwSetWindowSize(window_, w, h);
	}


	int Viewer::width() const {
		int w, h;
		glfwGetWindowSize(window_, &w, &h);
		return w;
	}


	int Viewer::height() const {
		int w, h;
		glfwGetWindowSize(window_, &w, &h);
		return h;
	}


	void Viewer::update() const {
		glfwPostEmptyEvent();
	}


	bool Viewer::mouse_press_event(int x, int y, int button, int modifiers) {
		if (modifiers == GLFW_MOD_SHIFT) {
			if (button == GLFW_MOUSE_BUTTON_LEFT) {
				bool found = false;
				const vec3& p = point_under_pixel(x, y, found);
				if (found) {
					camera_->setPivotPoint(p);
					// show, but hide the visual hint of pivot point after \p delay milliseconds.
					show_pivot_point_ = true;
					const int delay = 10000;
					Timer::single_shot(delay, [&]() {
						show_pivot_point_ = false;
						update();
					});

#if 1   // with animation
                    camera()->interpolateToLookAt(p);
#else   // without animation
                    const float coef = 0.1f;
                    const vec3& pos = coef * camera()->frame()->position() + (1.0f - coef)*p;
                    const quat& ori = camera()->frame()->orientation();
                    camera_->frame()->setPositionAndOrientation(pos, ori);
                    camera_->lookAt(p);
#endif
				}
				else {
					camera_->setPivotPoint(camera_->sceneCenter());
					show_pivot_point_ = false;
				}
			}
			else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
#if 1   // with animation
                camera()->interpolateToFitScene();
#else   // without animation
                camera()->showEntireScene();
#endif
                camera_->setPivotPoint(camera_->sceneCenter());
                show_pivot_point_ = false;
            }
        }

		camera_->frame()->action_start();
		return false;
	}


	bool Viewer::mouse_release_event(int x, int y, int button, int modifiers) {
		if (button == GLFW_MOUSE_BUTTON_LEFT && modifiers == GLFW_MOD_ALT) { // ZOOM_ON_REGION
			int xmin = std::min(mouse_pressed_x_, x);	int xmax = std::max(mouse_pressed_x_, x);
			int ymin = std::min(mouse_pressed_y_, y);	int ymax = std::max(mouse_pressed_y_, y);
			camera_->fitScreenRegion(xmin, ymin, xmax, ymax);
		}
		else
			camera_->frame()->action_end();

		button_ = -1;
		return false;
	}


	bool Viewer::mouse_drag_event(int x, int y, int dx, int dy, int button, int modifiers) {
		if (modifiers != GLFW_MOD_ALT) { // GLFW_MOD_ALT is reserved for zoom on region
			switch (button)
			{
			case GLFW_MOUSE_BUTTON_LEFT:
				camera_->frame()->action_rotate(x, y, dx, dy, camera_, modifiers == EASY3D_MOD_CONTROL);
				break;
			case GLFW_MOUSE_BUTTON_RIGHT:
				camera_->frame()->action_translate(x, y, dx, dy, camera_, modifiers == EASY3D_MOD_CONTROL);
				break;
			case GLFW_MOUSE_BUTTON_MIDDLE:
				if (dy != 0)
					camera_->frame()->action_zoom(dy > 0 ? 1 : -1, camera_);
				break;
			}
		}

		return false;
	}


	bool Viewer::mouse_free_move_event(int x, int y, int dx, int dy, int modifiers) {
		(void)x;
		(void)y;
		(void)dx;
		(void)dy;
		(void)modifiers;
		// highlight geometry primitives here
		return false;
	}


	bool Viewer::mouse_scroll_event(int x, int y, int dx, int dy) {
		(void)x;
		(void)y;
		(void)dx;
		camera_->frame()->action_zoom(dy, camera_);
		return false;
	}


	Model* Viewer::current_model() const {
		if (models_.empty())
			return nullptr;
		if (model_idx_ < models_.size())
			return models_[model_idx_];
		return nullptr;
	}

	bool Viewer::key_press_event(int key, int modifiers) {
		if (key == GLFW_KEY_F1 && modifiers == 0)
			LOG(INFO) << usage();

		else if (key == GLFW_KEY_LEFT && modifiers == 0) {
			float angle = static_cast<float>(1 * M_PI / 180.0); // turn left, 1 degrees each step
			camera_->frame()->action_turn(angle, camera_);
		}
		else if (key == GLFW_KEY_RIGHT && modifiers == 0) {
			float angle = static_cast<float>(1 * M_PI / 180.0); // turn right, 1 degrees each step
			camera_->frame()->action_turn(-angle, camera_);
		}

		else if (key == GLFW_KEY_UP && modifiers == 0) {	// move camera forward
			float step = 0.05f * camera_->sceneRadius();
			camera_->frame()->translate(camera_->frame()->inverseTransformOf(vec3(0.0, 0.0, -step)));
		}
		else if (key == GLFW_KEY_DOWN && modifiers == 0) {// move camera backward
			float step = 0.05f * camera_->sceneRadius();
			camera_->frame()->translate(camera_->frame()->inverseTransformOf(vec3(0.0, 0.0, step)));
		}

		else if (key == GLFW_KEY_LEFT && modifiers == EASY3D_MOD_CONTROL) {	// move camera left
			float step = 0.05f * camera_->sceneRadius();
			camera_->frame()->translate(camera_->frame()->inverseTransformOf(vec3(-step, 0.0, 0.0)));
		}
		else if (key == GLFW_KEY_RIGHT && modifiers == EASY3D_MOD_CONTROL) {	// move camera right
			float step = 0.05f * camera_->sceneRadius();
			camera_->frame()->translate(camera_->frame()->inverseTransformOf(vec3(step, 0.0, 0.0)));
		}

		else if (key == GLFW_KEY_UP && modifiers == EASY3D_MOD_CONTROL) {	// move camera up
			float step = 0.05f * camera_->sceneRadius();
			camera_->frame()->translate(camera_->frame()->inverseTransformOf(vec3(0.0, step, 0.0)));
		}
		else if (key == GLFW_KEY_DOWN && modifiers == EASY3D_MOD_CONTROL) {	// move camera down
			float step = 0.05f * camera_->sceneRadius();
			camera_->frame()->translate(camera_->frame()->inverseTransformOf(vec3(0.0, -step, 0.0)));
		}

		else if (key == GLFW_KEY_A && modifiers == 0) {
			if (drawable_axes_)
				drawable_axes_->set_visible(!drawable_axes_->is_visible());
		}
		else if (key == GLFW_KEY_C && modifiers == 0) {
			if (current_model())
				fit_screen(current_model());
		}
		else if (key == GLFW_KEY_F && modifiers == 0) {
			fit_screen();
		}

		else if (key == GLFW_KEY_M && modifiers == 0) {
#if 0
			// NOTE: switching on/off MSAA in this way will affect all viewers because OpenGL
			//       is a state machine. For multi-window applications, you have to call
			//		 glDisable()/glEnable() before the individual draw functions.
			if (samples_ > 0) {
				if (glIsEnabled(GL_MULTISAMPLE)) {
					glDisable(GL_MULTISAMPLE);
					LOG(INFO) << title_ + ": MSAA disabled";
				}
				else {
					glEnable(GL_MULTISAMPLE);
					LOG(INFO) << title_ + ": MSAA enabled";
				}
			}
#else
			if (dynamic_cast<SurfaceMesh*>(current_model())) {
				auto drawables = current_model()->triangles_drawables();
				for (auto d : drawables)
					d->set_smooth_shading(!d->smooth_shading());
			}
#endif
		}

		else if (key == GLFW_KEY_P && modifiers == 0) {
			if (camera_->type() == Camera::PERSPECTIVE)
				camera_->setType(Camera::ORTHOGRAPHIC);
			else
				camera_->setType(Camera::PERSPECTIVE);
		}
		else if (key == GLFW_KEY_SPACE && modifiers == 0) {
			// Aligns camera
			Frame frame;
			frame.setTranslation(camera_->pivotPoint());
			camera_->frame()->alignWithFrame(&frame, true);

			// Aligns frame
			//if (manipulatedFrame())
			//	manipulatedFrame()->alignWithFrame(camera_->frame());
		}

		else if (key == GLFW_KEY_O && modifiers == EASY3D_MOD_CONTROL)
			open();
		else if (key == GLFW_KEY_S && modifiers == EASY3D_MOD_CONTROL)
			save();

		else if (key == GLFW_KEY_MINUS && modifiers == EASY3D_MOD_CONTROL)
			camera_->frame()->action_zoom(-1, camera_);
		else if (key == GLFW_KEY_EQUAL && modifiers == EASY3D_MOD_CONTROL)
			camera_->frame()->action_zoom(1, camera_);

		if (key == GLFW_KEY_K && modifiers == GLFW_MOD_ALT) { // add key frame
			easy3d::Frame* frame = camera()->frame();
			camera()->keyFrameInterpolator()->addKeyFrame(*frame);
			// update scene bounding box to make sure the path is within the view frustum
			float old_radius = camera()->sceneRadius();
			float candidate_radius = distance(camera()->sceneCenter(), frame->position());
			camera()->setSceneRadius(std::max(old_radius, candidate_radius));
		}
		else if (key == GLFW_KEY_D && modifiers == EASY3D_MOD_CONTROL) { // delete path
			camera()->keyFrameInterpolator()->deletePath();

			// update scene bounding box
			Box3 box;
			for (auto m : models_)
				box.add_box(m->bounding_box());
			for (auto d : drawables_)
				box.add_box(d->bounding_box());
			camera_->setSceneBoundingBox(box.min(), box.max());
		}
		else if (key == GLFW_KEY_K && modifiers == EASY3D_MOD_CONTROL) { // play the path
			if (camera()->keyFrameInterpolator()->interpolationIsStarted())
				camera()->keyFrameInterpolator()->stopInterpolation();
			else
				camera()->keyFrameInterpolator()->startInterpolation();
		}
		else if (key == GLFW_KEY_T && modifiers == 0) {
			show_camera_path_ = !show_camera_path_;
		}

		else if (key == GLFW_KEY_LEFT_BRACKET && modifiers == 0) {
            for (auto m : models_) {
                for (auto d : m->lines_drawables()) {
                    float size = d->line_width() - 1.0f;
                    if (size < 1)
                        size = 1;
                    d->set_line_width(size);
                }
            }
        } if (key == GLFW_KEY_RIGHT_BRACKET && modifiers == 0) {
            for (auto m : models_) {
                for (auto d : m->lines_drawables()) {
                    float size = d->line_width() + 1.0f;
                    d->set_line_width(size);
                }
            }
        }

		else if (key == GLFW_KEY_MINUS && modifiers == 0) {
			for (auto m : models_) {
				for (auto d : m->points_drawables()) {
					float size = d->point_size() - 1.0f;
					if (size < 1)
						size = 1;
					d->set_point_size(size);
				}
			}
		}
		else if (key == GLFW_KEY_EQUAL && modifiers == 0) {
			for (auto m : models_) {
				for (auto d : m->points_drawables()) {
					float size = d->point_size() + 1.0f;
					d->set_point_size(size);
				}
			}
		}

		else if (key == GLFW_KEY_COMMA && modifiers == 0) {
			if (models_.empty())
				model_idx_ = -1;
			else
				model_idx_ = int((model_idx_ - 1 + models_.size()) % models_.size());
			if (model_idx_ >= 0) {
				fit_screen(models_[model_idx_]);
				LOG(INFO) << "current model: " << model_idx_ << ", " << models_[model_idx_]->name();
			}
		}
		else if (key == GLFW_KEY_PERIOD && modifiers == 0) {
			if (models_.empty())
				model_idx_ = -1;
			else
				model_idx_ = int((model_idx_ + 1) % models_.size());
			if (model_idx_ >= 0) {
				fit_screen(models_[model_idx_]);
				LOG(INFO) << "current model: " << model_idx_ << ", " << models_[model_idx_]->name();
			}
		}
		else if (key == GLFW_KEY_DELETE && modifiers == 0) {
			if (current_model())
				delete_model(current_model());
		}
		else if (key == GLFW_KEY_E && modifiers == 0) {
			if (current_model()) {
				LinesDrawable* drawable = current_model()->lines_drawable("edges");
				if (!drawable) {
					if (!dynamic_cast<PointCloud*>(current_model())) { // no default "edges" drawable for point clouds
						drawable = current_model()->add_lines_drawable("edges");
						renderer::update_buffer(current_model(), drawable);
					}
				}
				else
					drawable->set_visible(!drawable->is_visible());
			}
		}
		else if (key == GLFW_KEY_V && modifiers == 0) {
			if (current_model()) {
				PointsDrawable* drawable = current_model()->points_drawable("vertices");
				if (!drawable) {
					drawable = current_model()->add_points_drawable("vertices");
					renderer::update_buffer(current_model(), drawable);
                    drawable->set_impostor_type(PointsDrawable::SPHERE);
                    drawable->set_point_size(setting::surface_mesh_vertices_point_size);
				}
				else
					drawable->set_visible(!drawable->is_visible());
			}
		}

		else if (key == GLFW_KEY_B && modifiers == 0) {
			SurfaceMesh* mesh = dynamic_cast<SurfaceMesh*>(current_model());
			if (mesh) {
				auto drawable = mesh->lines_drawable("borders");
				if (!drawable) {
					auto prop = mesh->get_vertex_property<vec3>("v:point");
					std::vector<vec3> points;
					for (auto e : mesh->edges()) {
						if (mesh->is_boundary(e)) {
							points.push_back(prop[mesh->vertex(e, 0)]);
							points.push_back(prop[mesh->vertex(e, 1)]);
						}
					}
					if (!points.empty()) {
						drawable = mesh->add_lines_drawable("borders");
						drawable->update_vertex_buffer(points);
						drawable->set_default_color(setting::surface_mesh_borders_color);
						drawable->set_per_vertex_color(false);
						drawable->set_impostor_type(LinesDrawable::CYLINDER);
						drawable->set_line_width(setting::surface_mesh_borders_line_width);
					}
				}
				else
					drawable->set_visible(!drawable->is_visible());
			}
		}

        else if (key == GLFW_KEY_L && modifiers == 0) { // locked vertices
            SurfaceMesh* mesh = dynamic_cast<SurfaceMesh*>(current_model());
            if (mesh) {
                auto drawable = mesh->points_drawable("locks");
                if (!drawable) {
                    auto lock = mesh->get_vertex_property<bool>("v:lock");
                    if (lock) {
                        auto prop = mesh->get_vertex_property<vec3>("v:point");
                        std::vector<vec3> points;
                        for (auto v : mesh->vertices()) {
                            if (lock[v])
                                points.push_back(prop[v]);
                        }
                        if (!points.empty()) {
                            drawable = mesh->add_points_drawable("locks");
                            drawable->update_vertex_buffer(points);
                            drawable->set_default_color(vec4(1, 1, 0, 1));
                            drawable->set_per_vertex_color(false);
                            drawable->set_impostor_type(PointsDrawable::SPHERE);
                            drawable->set_point_size(setting::surface_mesh_vertices_point_size + 3);
                        }
                    }
                }
                else
                    drawable->set_visible(!drawable->is_visible());
            }
        }

		else if (key == GLFW_KEY_D && modifiers == 0) {
			if (current_model()) {
				std::cout << "----------- " << file_system::simple_name(current_model()->name()) << " -----------\n";
				if (dynamic_cast<SurfaceMesh*>(current_model())) {
					auto model = dynamic_cast<SurfaceMesh*>(current_model());
					std::cout << "model is a surface mesh. #face: " << std::to_string(model->n_faces())
							  << ", #vertex: " + std::to_string(model->n_vertices())
							  << ", #edge: " + std::to_string(model->n_edges()) << std::endl;
				}
				else if (dynamic_cast<PointCloud*>(current_model())) {
					auto model = dynamic_cast<PointCloud*>(current_model());
					std::cout << "model is a point cloud. #vertex: " + std::to_string(model->n_vertices()) << std::endl;
				}
				else if (dynamic_cast<Graph*>(current_model())) {
					auto model = dynamic_cast<Graph*>(current_model());
					std::cout << "model is a graph. #vertex: " + std::to_string(model->n_vertices())
							  << ", #edge: " + std::to_string(model->n_edges()) << std::endl;
				}
				if (!current_model()->points_drawables().empty()) {
					std::cout << "points drawables:\n";
					for (auto d : current_model()->points_drawables())
						d->drawable_stats();
				}
				if (!current_model()->lines_drawables().empty()) {
					std::cout << "lines drawables:\n";
					for (auto d : current_model()->lines_drawables())
                        d->drawable_stats();
				}
				if (!current_model()->triangles_drawables().empty()) {
					std::cout << "triangles drawables:\n";
					for (auto d : current_model()->triangles_drawables())
                        d->drawable_stats();
				}

				current_model()->property_stats();
			}
		}

		else if (key == GLFW_KEY_R && modifiers == 0) {
			// Reload the shader(s) - useful for writing/debugging shader code.
			ShaderManager::reload();
		}

		else if (key == GLFW_KEY_S && modifiers == 0) {
			snapshot();
		}

		else if (key == GLFW_KEY_F4 && modifiers == GLFW_MOD_ALT) {
			glfwSetWindowShouldClose(window_, true);
		}

		return false;
	}


	bool Viewer::key_release_event(int key, int modifiers) {
		(void)key;
		(void)modifiers;
		return false;
	}


	bool Viewer::char_input_event(unsigned int codepoint) {
		(void)codepoint;
		//switch (codepoint) {
		//case '-':	break;
		//case 'c':	break;
		//case 'C':	break;
		//case '0':	break;
		//default:
		//	return false;
		//}

		return false;
	}


	bool Viewer::drop_event(const std::vector<std::string> & filenames) {
		int count = 0;
		for (auto& name : filenames) {
			if (add_model(name))
				++count;
		}

		if (count > 0) {
			model_idx_ = static_cast<int>(models_.size()) - 1; // make the last one current
			fit_screen();
			return true;
		}
		return false;
	}


	vec3 Viewer::point_under_pixel(int x, int y, bool &found) const {
		// GLFW (same as Qt) uses upper corner for its origin while GL uses the lower corner.
		int glx = x;
		int gly = height() - 1 - y;

		// NOTE: when dealing with OpenGL, we alway work in the highdpi screen space
#if defined(__APPLE__)
		glx = static_cast<int>(glx * dpi_scaling());
		gly = static_cast<int>(gly * dpi_scaling());
#endif
		float depth = std::numeric_limits<float>::max();
		glPixelStorei(GL_PACK_ALIGNMENT, 1);				easy3d_debug_log_gl_error;
		glReadPixels(glx, gly, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);					easy3d_debug_log_gl_error;
		found = depth < 1.0f;
		if (found) {
			vec3 point(float(x), float(y), depth);
			// The input to unprojectedCoordinatesOf() is defined in the screen coordinate system
			point = camera_->unprojectedCoordinatesOf(point);
			return point;
		}
		else
			return vec3(0);
	}


	inline double get_seconds() {
		return std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
	}


	void Viewer::run() {
        // make sure scene fits the screen when the window appears
        fit_screen();

        // initialize before showing the window because it can be slow
		init();

		// show the window
		glfwShowWindow(window_);

		// TODO: make it member variable
		bool is_animating = false;

		try {
			// Rendering loop
			const int num_extra_frames = 5;
			const double animation_max_fps = 30;
			int frame_counter = 0;

			while (!glfwWindowShouldClose(window_)) {
				if (!glfwGetWindowAttrib(window_, GLFW_VISIBLE)) // not visible
					continue;

				double tic = get_seconds();
				pre_draw();

				gpu_timer_->start();
				draw();
				gpu_timer_->stop();
				gpu_time_ = gpu_timer_->time();

				post_draw();
				glfwSwapBuffers(window_);

				if (is_animating || frame_counter++ < num_extra_frames)
				{
					glfwPollEvents();
					// In microseconds
					double duration = 1000000.* (get_seconds() - tic);
					const double min_duration = 1000000. / animation_max_fps;
					if (duration < min_duration)
						std::this_thread::sleep_for(std::chrono::microseconds(static_cast<int>(min_duration - duration)));
				}
				else {
					/* Wait for mouse/keyboard or empty refresh events */
					glfwWaitEvents();
					frame_counter = 0;
				}
			}

			/* Process events once more */
			glfwPollEvents();
		}
		catch (const std::exception &e) {
			LOG(ERROR) << "Caught exception in main loop: " << e.what();
		}

		cleanup();
	}


	void Viewer::init() {
		std::cerr << usage() << std::endl;

		glEnable(GL_DEPTH_TEST);

		glClearDepthf(1.0f);
		glClearColor(background_color_[0], background_color_[1], background_color_[2], background_color_[3]);
	}


	std::string Viewer::usage() const {
		return std::string(
                " ------------------------------------------------------------------\n"
                " Easy3D viewer usage:                                              \n"
                " ------------------------------------------------------------------\n"
                "  F1:                  Help                                        \n"
                " ------------------------------------------------------------------\n"
                "  Ctrl + 'o':          Open file                                   \n"
                "  Ctrl + 's':          Save file                                   \n"
                "  Fn + Delete:         Delete current model                        \n"
                "  '<' or '>':          Switch between models                       \n"
                "  's':                 Snapshot                                    \n"
                " ------------------------------------------------------------------\n"
                "  'p':                 Toggle perspective/orthographic projection)	\n"
                "  Left mouse:          Orbit-rotate the camera                     \n"
                "  Right mouse:         Pan-move the camera                         \n"
                "  Ctrl + Left mouse:   Orbit-rotate the camera (screen based)      \n"
                "  Ctrl + Right mouse:  Pan-move camera vertically or horizontally  \n"
                "  Middle or Wheel:     Zoom in/out                                 \n"
                "  Ctrl + '+'/'-':      Zoom in/out                                 \n"
                "  Left/Right           Turn camera left/right                      \n"
                "  Ctrl + Left/Right:   Move camera left/right                      \n"
                "  Up/Down:             Move camera forward/backward                \n"
                "  Ctrl + Up/Down:      Move camera up/down                         \n"
                " ------------------------------------------------------------------\n"
                "  'f':                 Fit screen (all models)                     \n"
                "  'c':                 Fit screen (current model only)             \n"
                "  Shift + Left mouse:  Zoom to target point on model               \n"
                "  Shift + Right mouse: Zoom too fit screen                         \n"
                "  Alt + Left mouse:    Zoom too fit region                         \n"
                " ------------------------------------------------------------------\n"
                "  '-'/'=':             Decrease/Increase point size                \n"
                "  '{'/'}':             Decrease/Increase line width                \n"
                "  'a':                 Toggle axes									\n"
                "  'b':                 Toggle borders								\n"
                "  'e':                 Toggle edges							    \n"
                "  'v':                 Toggle vertices                             \n"
                "  'm':                 Toggle smooth shading (for SurfaceMesh)     \n"
                "  'd':                 Print model info (drawables, properties)    \n"
                " ------------------------------------------------------------------\n"
		);
	}


	Model* Viewer::add_model(const std::string &file_name, bool create_default_drawables) {
		for (auto m : models_) {
			if (m->name() == file_name) {
				LOG(WARNING) << "model has already been added to the viewer: " << file_name;
				return m;
			}
		}

		const std::string& ext = file_system::extension(file_name, true);
		Model* model = nullptr;
		if (ext == "obj") { // mesh
			model = SurfaceMeshIO::load(file_name);
		}

		if (model) {
			model->set_name(file_name);
			add_model(model, create_default_drawables);
		}
		return model;
	}


	void Viewer::create_drawables(Model *model) {
		if (dynamic_cast<PointCloud*>(model)) {
			PointCloud* cloud = dynamic_cast<PointCloud*>(model);
            renderer::update_buffer(cloud, cloud->add_points_drawable("vertices"));
		}
		else if (dynamic_cast<SurfaceMesh*>(model)) {
			SurfaceMesh* mesh = dynamic_cast<SurfaceMesh*>(model);
            renderer::update_buffer(mesh, mesh->add_triangles_drawable("faces"));

            if (setting::surface_mesh_show_edges)
                renderer::update_buffer(mesh, mesh->add_lines_drawable("edges"));

            if (setting::surface_mesh_show_vertices)
                renderer::update_buffer(mesh, mesh->add_points_drawable("edges"));
		}
		else if (dynamic_cast<Graph*>(model)) {
			Graph* graph = dynamic_cast<Graph*>(model);

			// create points drawable for the edges
			PointsDrawable* vertices = graph->add_points_drawable("vertices");
			renderer::update_buffer(graph, vertices);
            vertices->set_default_color(setting::graph_vertices_color);
            vertices->set_point_size(setting::graph_vertices_point_size);
            vertices->set_impostor_type(PointsDrawable::SPHERE);

			// create liens drawable for the edges
			LinesDrawable* edges = graph->add_lines_drawable("edges");
			renderer::update_buffer(graph, edges);
            edges->set_default_color(setting::graph_edges_color);
            edges->set_line_width(setting::graph_edges_line_width);
            edges->set_impostor_type(LinesDrawable::CYLINDER);
		}
	}


	Model* Viewer::add_model(Model *model, bool create_default_drawables) {
		if (!model) {
			LOG(WARNING) << "model is NULL.";
			return nullptr;
		}
		for (auto m : models_) {
			if (model == m) {
				LOG(WARNING) << "model has already been added to the viewer: " << m->name();
				return nullptr;
			}
		}

		if (create_default_drawables) {
			if (model->n_vertices() > 0) {
                create_drawables(model);
            }
			else
				LOG(WARNING) << "drawable cannot be created due to no vertices.";
		}

		int pre_idx = model_idx_;
		models_.push_back(model);
		model_idx_ = static_cast<int>(models_.size()) - 1; // make the last one current

		return model;
	}


	bool Viewer::delete_model(Model* model) {
		if (!model) {
			LOG(WARNING) << "model is NULL.";
			return false;
		}

		auto pos = std::find(models_.begin(), models_.end(), model);
		if (pos != models_.end()) {
			int pre_idx = model_idx_;
			const std::string name = model->name();
			models_.erase(pos);
			delete model;
			model_idx_ = static_cast<int>(models_.size()) - 1; // make the last one current
			LOG(INFO) << "model deleted: " << name;

			if (model_idx_ != pre_idx) {
				if (model_idx_ >= 0)
					LOG(INFO) << "current model: " << model_idx_ << ", " << models_[model_idx_]->name();
			}
			return true;
		}
		else {
			LOG(WARNING) << "no such model: " << model->name();
			return false;
		}
	}


	bool Viewer::add_drawable(Drawable* drawable) {
		if (!drawable) {
			LOG(WARNING) << "drawable is NULL.";
			return false;
		}
		for (auto d : drawables_) {
			if (drawable == d) {
				LOG(WARNING) << "drawable has already been added to the viewer.";
				return false;
			}
		}

		drawables_.push_back(drawable);
		return true;
	}


	bool Viewer::delete_drawable(Drawable* drawable) {
		if (!drawable) {
			LOG(WARNING) << "drawable is NULL.";
			return false;
		}

		auto pos = std::find(drawables_.begin(), drawables_.end(), drawable);
		if (pos != drawables_.end()) {
			drawables_.erase(pos);
			delete drawable;
			return true;
		}
		else {
			LOG(WARNING) << "no such drawable: " << drawable->name();
			return false;
		}
	}


	void Viewer::fit_screen(const easy3d::Model *model) {
		if (!model && models_.empty() && drawables_.empty())
			return;

		auto visual_box = [](const Model* m) -> Box3 {
			Box3 box = m->bounding_box();
			for (auto d : m->points_drawables())    box.add_box(d->bounding_box());
			for (auto d : m->lines_drawables())     box.add_box(d->bounding_box());
			for (auto d : m->triangles_drawables()) box.add_box(d->bounding_box());
			return box;
		};

		Box3 box;
		if (model)
			box = visual_box(model);
		else {
			for (auto m : models_)
				box.add_box(visual_box(m));
			for (auto d : drawables_)
				box.add_box(d->bounding_box());
		}

		if (box.initialized()) {
			camera_->setSceneBoundingBox(box.min(), box.max());
			camera_->showEntireScene();
			update();
		}
	}


	bool Viewer::open() {
		const std::string title("Please choose a file");
		const std::string& default_path = resource::directory() + "/data/";
		const std::vector<std::string>& filters = {
			"Mesh Files (*.obj *.ply *.off *.stl *.poly)" , "*.obj *.ply *.off *.stl *.poly" ,
			"Point Cloud Files (*.bin *.ply *.xyz *.bxyz *.las *.laz *.vg *.bvg *.ptx)", "*.bin *.ply *.xyz *.bxyz *.las *.laz *.vg *.bvg *.ptx",
			"All Files (*.*)", "*"
		};
		const std::vector<std::string>& file_names = dialog::open(title, default_path, filters, true);

		int count = 0;
		for (const auto& file_name : file_names) {
			if (add_model(file_name))
				++count;
		}

		if (count > 0) {
			fit_screen();
			return true;
		}
		return false;
	}


	bool Viewer::save() const {
		const Model* m = current_model();
		if (!m) {
			LOG(ERROR) << "no model exists";
			return false;
		}

		const std::string& title = "Please choose a file name";
        const std::vector<std::string> &filters = {
                "Point Cloud Files (*.xyz)", "*.xyz",
                "Mesh Files (*.obj)", "*.obj",
                "All Files (*.*)", "*"
        };

		std::string default_file_name = m->name();
		const bool warn_overwrite = true;
		const std::string& file_name = dialog::save(title, default_file_name, filters, warn_overwrite);
		if (file_name.empty())
			return false;

        const std::string& ext = file_system::extension(file_name, true);

		bool saved = false;
		if (dynamic_cast<const SurfaceMesh*>(m)) {
            saved = SurfaceMeshIO::save(file_name, dynamic_cast<const SurfaceMesh *>(m));
            if (ext != "obj") {
                LOG(ERROR) << "surface mesh file must have extension 'obj'";
                return false;
            }
        }
        else if (dynamic_cast<const PointCloud*>(m)) {
            if (ext != "xyz") {
                LOG(ERROR) << "point cloud file must have extension 'xyz'";
                return false;
            }
            saved = save_xyz(file_name, *dynamic_cast<const PointCloud *>(m));
        }
		if (saved)
			LOG(INFO) << "file successfully saved";

		return saved;
	}


	bool Viewer::snapshot(bool bk_white /* = true*/) const {
		const std::string& title = "Please choose a file name";
		std::string default_file_name("untitled.png");
		if (current_model())
			default_file_name = file_system::replace_extension(current_model()->name(), "png");
		const std::vector<std::string>& filters = {
			"Image Files (*.png)" , "*.png",
			"All Files (*.*)", "*"
		};

		const bool warn_overwrite = true;
		const std::string& file_name = dialog::save(title, default_file_name, filters, warn_overwrite);
		if (file_name.empty())
			return false;

		const std::string& ext = file_system::extension(file_name, true);
		if (ext != "png") {
			LOG(ERROR) << "snapshot format must be 'png'";
			return false;
		}

		int w, h;
		glfwGetFramebufferSize(window_, &w, &h);
		FramebufferObject fbo(w, h, samples_);
		fbo.add_color_buffer();
		fbo.add_depth_buffer();

		fbo.bind();

		if (bk_white)
			glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		else
			glClearColor(background_color_[0], background_color_[1], background_color_[2], background_color_[3]);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

		const_cast<Viewer*>(this)->draw();

		fbo.release();

#if 1
		// also saves the camera parameters
		std::ofstream output(file_name + ".camera");
		output << "position:      " << camera_->position() << std::endl;
        output << "orientation:   " << camera_->orientation() << std::endl;
        output << "field_of_view: " << camera_->fieldOfView() << std::endl;
#endif

#if 1   // color render buffer
		return fbo.snapshot_color(0, file_name);
#else
		// depth buffer
		return fbo.snapshot_depth(file_name);
#endif
	}


	void Viewer::draw_corner_axes() {
		ShaderProgram* program = ShaderManager::get_program("surface/surface_color");
		if (!program) {
			std::vector<ShaderProgram::Attribute> attributes;
			attributes.emplace_back(ShaderProgram::Attribute(ShaderProgram::POSITION, "vtx_position"));
			attributes.emplace_back(ShaderProgram::Attribute(ShaderProgram::COLOR, "vtx_color"));
			attributes.emplace_back(ShaderProgram::Attribute(ShaderProgram::NORMAL, "vtx_normal"));
			program = ShaderManager::create_program_from_files("surface/surface_color", attributes);
		}
		if (!program)
			return;

		if (!drawable_axes_) {
			float base = 0.5f;   // the cylinder length, relative to the allowed region
			float head = 0.2f;   // the cone length, relative to the allowed region
			std::vector<vec3> points, normals, colors;
			opengl::prepare_cylinder(0.03, 10, vec3(0, 0, 0), vec3(base, 0, 0), vec3(1, 0, 0), points, normals, colors);
			opengl::prepare_cylinder(0.03, 10, vec3(0, 0, 0), vec3(0, base, 0), vec3(0, 1, 0), points, normals, colors);
			opengl::prepare_cylinder(0.03, 10, vec3(0, 0, 0), vec3(0, 0, base), vec3(0, 0, 1), points, normals, colors);
			opengl::prepare_cone(0.06, 20, vec3(base, 0, 0), vec3(base + head, 0, 0), vec3(1, 0, 0), points, normals, colors);
			opengl::prepare_cone(0.06, 20, vec3(0, base, 0), vec3(0, base + head, 0), vec3(0, 1, 0), points, normals, colors);
			opengl::prepare_cone(0.06, 20, vec3(0, 0, base), vec3(0, 0, base + head), vec3(0, 0, 1), points, normals, colors);
			opengl::prepare_sphere(vec3(0, 0, 0), 0.06, 20, 20, vec3(0, 1, 1), points, normals, colors);
			drawable_axes_ = new TrianglesDrawable("corner_axes");
			drawable_axes_->update_vertex_buffer(points);
			drawable_axes_->update_normal_buffer(normals);
			drawable_axes_->update_color_buffer(colors);
			drawable_axes_->set_per_vertex_color(true);
		}
		if (!drawable_axes_->is_visible())
			return;

		// The viewport and the scissor are changed to fit the lower left corner.
		int viewport[4], scissor[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		glGetIntegerv(GL_SCISSOR_BOX, scissor);

		static int corner_frame_size = static_cast<int>(100 * dpi_scaling());
		glViewport(0, 0, corner_frame_size, corner_frame_size);
		glScissor(0, 0, corner_frame_size, corner_frame_size);

		// To make the axis appear over other objects: reserve a tiny bit of the
		// front depth range. NOTE: do remember to restore it later.
		glDepthRangef(0, 0.001f);

		const mat4& proj = transform::ortho(-1, 1, -1, 1, -1, 1);
		const mat4& view = camera_->orientation().inverse().matrix();
        const mat4 &MVP = proj * view;

        // camera position is defined in world coordinate system.
        const vec3 &wCamPos = camera_->position();
        // it can also be computed as follows:
        //const vec3& wCamPos = invMV * vec4(0, 0, 0, 1);
        const mat4 &MV = camera_->modelViewMatrix();
        const vec4 &wLightPos = inverse(MV) * setting::light_position;

        program->bind();
        program->set_uniform("MVP", MVP)
                ->set_uniform("lighting", true)
                ->set_uniform("two_sides_lighting", false)
                ->set_uniform("smooth_shading", true)
                ->set_uniform("wLightPos", wLightPos)
                ->set_uniform("wCamPos", wCamPos)
                ->set_uniform("ssaoEnabled", false)
                ->set_uniform("per_vertex_color", true)
                ->set_uniform("distinct_back_color", false)
                ->set_block_uniform("Material", "ambient", setting::material_ambient)
                ->set_block_uniform("Material", "specular", setting::material_specular)
                ->set_block_uniform("Material", "shininess", &setting::material_shininess)
                ->set_uniform("hightlight_id_min", -1)
                ->set_uniform("hightlight_id_max", -1);
        drawable_axes_->gl_draw(false);
        program->release();

        // restore
        glScissor(scissor[0], scissor[1], scissor[2], scissor[3]);
        glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
        glDepthRangef(0.0f, 1.0f);
    }


	void Viewer::pre_draw() {
		glfwMakeContextCurrent(window_);
		glClearColor(background_color_[0], background_color_[1], background_color_[2], 1.0f);
		glClearDepth(1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	}


	void Viewer::post_draw() {
		// shown only when it is not animating
		if (show_camera_path_ && !camera()->keyFrameInterpolator()->interpolationIsStarted())
			camera()->draw_paths();

		if (show_pivot_point_) {
			ShaderProgram* program = ShaderManager::get_program("lines/lines_plain_color");
			if (!program) {
				std::vector<ShaderProgram::Attribute> attributes;
				attributes.emplace_back(ShaderProgram::Attribute(ShaderProgram::POSITION, "vtx_position"));
				attributes.emplace_back(ShaderProgram::Attribute(ShaderProgram::COLOR, "vtx_color"));
				program = ShaderManager::create_program_from_files("lines/lines_plain_color", attributes);
			}
			if (!program)
				return;

#if defined(__APPLE__)
			const float size = 10;
#else
			const float size = static_cast<float>(10 * dpi_scaling());
#endif
			LinesDrawable drawable("pivot_point");
			const vec3& pivot = camera()->projectedCoordinatesOf(camera()->pivotPoint());
			std::vector<vec3> points = {
					vec3(pivot.x - size, pivot.y, 0.5f), vec3(pivot.x + size, pivot.y, 0.5f),
					vec3(pivot.x, pivot.y - size, 0.5f), vec3(pivot.x, pivot.y + size, 0.5f)
			};
			drawable.update_vertex_buffer(points);

			const mat4& proj = transform::ortho(0.0f, static_cast<float>(width()), static_cast<float>(height()), 0.0f, 0.0f, -1.0f);
			glDisable(GL_DEPTH_TEST);   // always on top
			program->bind();
			program->set_uniform("MVP", proj);
			program->set_uniform("per_vertex_color", false);
			program->set_uniform("default_color", vec3(0.0f, 0.0f, 1.0f));
			drawable.gl_draw(false);
			program->release();
			glEnable(GL_DEPTH_TEST);   // restore
		}

		if (button_ == GLFW_MOUSE_BUTTON_LEFT && modifiers_ == GLFW_MOD_ALT) {
		    const Rect rect(mouse_pressed_x_, mouse_x_, mouse_pressed_y_, mouse_y_);
            if (rect.width() > 0 || rect.height() > 0) {
                // draw the boundary of the rect
                opengl::draw_quad_wire(rect, vec4(0.0f, 0.0f, 1.0f, 1.0f), width(), height(), -1.0f);
                // draw the transparent face
                glEnable(GL_BLEND);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                opengl::draw_quad_filled(rect, vec4(0.0f, 0.0f, 1.0f, 0.2f), width(), height(), -0.9f);
                glDisable(GL_BLEND);
            }
        }

		draw_corner_axes();
	}


	void Viewer::draw() const {
		for (const auto m : models_) {
			if (!m->is_visible())
				continue;

			for (auto d : m->points_drawables()) {
				if (d->is_visible())
					d->draw(camera(), false);
			}

			// Let's check if edges and surfaces are both shown. If true, we
			// make the depth coordinates of the surface smaller, so that displaying
			// the mesh and the surface together does not cause Z-fighting.
			std::size_t count = 0;
			for (auto d : m->lines_drawables()) {
				if (d->is_visible()) {
					d->draw(camera(), false);
					++count;
				}
			}

			if (count > 0) {
				glEnable(GL_POLYGON_OFFSET_FILL);
				glPolygonOffset(0.5f, -0.0001f);
			}
			for (auto d : m->triangles_drawables()) {
				if (d->is_visible())
					d->draw(camera(), false);
			}
			if (count > 0)
				glDisable(GL_POLYGON_OFFSET_FILL);

		}

		for (auto d : drawables_) {
			if (d->is_visible())
				d->draw(camera(), false);
		}
	}


}
