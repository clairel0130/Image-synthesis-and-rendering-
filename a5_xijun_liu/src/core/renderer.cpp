/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/core.h>
#include <core/accel.h>
#include <core/renderer.h>
#include <GL/glew.h>

#ifdef __APPLE__
#include "SDL.h"
#include <OpenGL/gl.h>
#else
#ifdef _WIN32
#include <GL/gl.h>
#include "SDL.h"
#else
#include <GL/gl.h>
#include "SDL2/SDL.h"
#endif
#endif

#include <bsdfs/diffuse.h>

#include <integrators/normal.h>
#include <renderpasses/normal.h>

#include <integrators/simple.h>
#include <renderpasses/simple.h>
#include <bsdfs/phong.h>

#include <integrators/ao.h>
#include <integrators/ro.h>
#include <renderpasses/ssao.h>

#include <integrators/direct.h>

#include <integrators/path.h>
#include <renderpasses/gi.h>
#include <bsdfs/mixture.h>


TR_NAMESPACE_BEGIN

Renderer::Renderer(const Config& config) : scene(config) { }

bool Renderer::init(const bool isRealTime, bool nogui) {
    realTime = isRealTime;
    this->nogui = nogui;
    realTimeCameraFree = false;

    if (!scene.load(isRealTime)) return false;

    if (realTime) {
        if (scene.config.renderpass == ENormalRenderPass) {
            renderpass = std::unique_ptr<NormalPass>(new NormalPass(scene));
        }
        else if (scene.config.renderpass == EDirectRenderPass) {
            renderpass = std::unique_ptr<SimplePass>(new SimplePass(scene));
        }
        else if (scene.config.renderpass == ESSAORenderPass) {
            renderpass = std::unique_ptr<SSAOPass>(new SSAOPass(scene));
        }
        else if (scene.config.renderpass == EGIRenderPass) {
            renderpass = std::unique_ptr<GIPass>(new GIPass(scene));
        }
        else {
            throw std::runtime_error("Invalid renderpass type");
        }

        bool succ = renderpass.get()->initOpenGL(scene.config.width, scene.config.height);
        if (!succ) return false;

        return renderpass->init(scene.config);
    } else {
        if (scene.config.integrator == ENormalIntegrator) {
            integrator = std::unique_ptr<NormalIntegrator>(new NormalIntegrator(scene));
        }
        else if (scene.config.integrator == EAOIntegrator) {
            integrator = std::unique_ptr<AOIntegrator>(new AOIntegrator(scene));
        } else if (scene.config.integrator == EROIntegrator) {
            integrator = std::unique_ptr<ROIntegrator>(new ROIntegrator(scene));
        }
        else if (scene.config.integrator == ESimpleIntegrator) {
            integrator = std::unique_ptr<SimpleIntegrator>(new SimpleIntegrator(scene));
        }
        else if (scene.config.integrator == EDirectIntegrator) {
            integrator = std::unique_ptr<DirectIntegrator>(new DirectIntegrator(scene));
        }
        else if (scene.config.integrator == EPathTracerIntegrator) {
            integrator = std::unique_ptr<PathTracerIntegrator>(new PathTracerIntegrator(scene));
        }
        else {
            throw std::runtime_error("Invalid integrator type");
        }

        return integrator->init();
    }
}

void Renderer::render() {
    if (realTime) {
        // TODO: Add previous assignment code (if needed)
        while(1) {
            SDL_Event event;
            SDL_PollEvent(&event);
            if (event.type == SDL_QUIT) {
                break;
            }
            renderpass->render();
            SDL_GL_SwapWindow(renderpass->window);
        }
    } else {
        /**
         * 1) Calculate the camera perspective, the camera-to-world transformation matrix and the aspect ratio.
         * 2) Clear integral RGB buffer.
         * 3) Loop over all pixels on the image plane.
         * 4) Generate a ray through each pixel center.
         * 5) Splat their contribution onto the image plane.
         */
        // TODO: Add previous assignment code (if needed)

        // TODO: Implement this
        // Camera perspective
        glm::mat4 CameraMatrix = glm::lookAt(scene.config.camera.o, scene.config.camera.at, scene.config.camera.up);
        float rad = (scene.config.camera.fov * 2.f * M_PI) / 360.f;
        float scaling = tan(rad / 2.f);
        float aspectRatio = (float)scene.config.width / (float) scene.config.height;
        float px, py;

        int i = 0;
        Sampler sample = TinyRender::Sampler(260654285);
        integrator->rgb->clear();

//        ThreadPool::ParallelFor(0, scene.config.height, [&] (int y) {
//            for (int x = 0; x < scene.config.width; ++x) {
//                v3f pixelColor = v3f(0.f, 0.f, 0.f);
//                for(size_t n = 0; n < scene.config.spp; n++){
//                    float r;
//                    if (scene.config.spp == 1){
//                        r = 0.5;
//                    }else{
//                        r = sample.next();
//                    }
//                    px = (x - scene.config.width / 2.f + r) / (scene.config.width / 2.f) * aspectRatio * scaling;
//                    py = -((y - scene.config.height / 2.f + r) / (scene.config.height / 2.f) * scaling);
//                    v4f AugmentM = v4f(px, py, -1.f, 0.f);
//                    v4f dir = AugmentM * CameraMatrix;
//                    Ray ray = TinyRender::Ray(scene.config.camera.o, dir);
//                    pixelColor = pixelColor + integrator->render(ray, sample);
//                }
//                pixelColor = pixelColor/scene.config.spp;
//                integrator->rgb->data[y * scene.config.width + x] = pixelColor;
//            }
//        });

        for (size_t y = 0; y < scene.config.height; y++) {
            for (size_t x = 0; x < scene.config.width; x++) {
                v3f pixelColor = v3f(0.f, 0.f, 0.f);
                for(size_t n = 0; n < scene.config.spp; n++){
                    float r;
                    if (scene.config.spp == 1){
                        r = 0.5;
                    }else{
                        r = sample.next();
                    }
                    px = (x - scene.config.width / 2.f + r) / (scene.config.width / 2.f) * aspectRatio * scaling;
                    py = -((y - scene.config.height / 2.f + r) / (scene.config.height / 2.f) * scaling);
                    v4f AugmentM = v4f(px, py, -1.f, 0.f);
                    v4f dir = AugmentM * CameraMatrix;
                    Ray ray = TinyRender::Ray(scene.config.camera.o, dir);
                    pixelColor = pixelColor + integrator->render(ray, sample);
                }
                pixelColor = pixelColor/scene.config.spp;
                integrator->rgb->data[i] = pixelColor;
                i++;
            }
        }
    }
}
/**
 * Post-rendering step.
 */
void Renderer::cleanUp() {
    if (realTime) {
        renderpass->cleanUp();
    } else {
        integrator->cleanUp();
    }
}

BSDF::BSDF(const WorldData& d, const Config& c, const size_t matID) : worldData(d), config(c) {
    emission = glm::make_vec3(worldData.materials[matID].emission);
}

Scene::Scene(const Config& config) : config(config) { }

bool Scene::load(bool isRealTime) {
    fs::path file(config.objFile);
    bool ret = false;
    std::string err;

    if (!file.is_absolute())
        file = (config.tomlFile.parent_path() / file).make_preferred();

    tinyobj::attrib_t* attrib_ = &worldData.attrib;
    std::vector<tinyobj::shape_t>* shapes_ = &worldData.shapes;
    std::vector<tinyobj::material_t>* materials_ = &worldData.materials;
    std::string* err_ = &err;
    const string filename_ = file.string();
    const string mtl_basedir_ = file.make_preferred().parent_path().string();
    ret = tinyobj::LoadObj(attrib_, shapes_, materials_, err_, filename_.c_str(), mtl_basedir_.c_str(), true);

    if (!err.empty()) { std::cout << "Error: " << err.c_str() << std::endl; }
    if (!ret) {
        std::cout << "Failed to load scene " << config.objFile << " " << std::endl;
        return false;
    }

    // Build list of BSDFs
    bsdfs = std::vector<std::unique_ptr<BSDF>>(worldData.materials.size());
    for (size_t i = 0; i < worldData.materials.size(); i++) {
        if (worldData.materials[i].illum == 7)
            bsdfs[i] = std::unique_ptr<BSDF>(new DiffuseBSDF(worldData, config, i));
        if (worldData.materials[i].illum != 5 && worldData.materials[i].illum != 7 && worldData.materials[i].illum != 8)
            bsdfs[i] = std::unique_ptr<BSDF>(new PhongBSDF(worldData, config, i));
        if (worldData.materials[i].illum == 8)
            bsdfs[i] = std::unique_ptr<BSDF>(new MixtureBSDF(worldData, config, i));
    }

    // Build list of emitters (and print what has been loaded)
    std::string nbShapes = worldData.shapes.size() > 1 ? " shapes" : " shape";
    std::cout << "Found " << worldData.shapes.size() << nbShapes << std::endl;
    worldData.shapesCenter.resize(worldData.shapes.size());
    worldData.shapesAABOX.resize(worldData.shapes.size());

    for (size_t i = 0; i < worldData.shapes.size(); i++) {
        const tinyobj::shape_t& shape = worldData.shapes[i];
        const BSDF* bsdf = bsdfs[shape.mesh.material_ids[0]].get();
        std::cout << "Mesh " << i << ": " << shape.name << " ["
                  << shape.mesh.indices.size() / 3 << " primitives | ";

        if (bsdf->isEmissive()) {
            Distribution1D faceAreaDistribution;
            float shapeArea = getShapeArea(i, faceAreaDistribution);
            emitters.emplace_back(Emitter{i, shapeArea, bsdf->emission, faceAreaDistribution});
            std::cout << "Emitter]" << std::endl;
        } else {
            std::cout << bsdf->toString() << "]" << std::endl;
        }

        // Build world AABB and shape centers
        worldData.shapesCenter[i] = v3f(0.0);
        for (auto idx: shape.mesh.indices) {
            v3f p = {worldData.attrib.vertices[3 * idx.vertex_index + 0],
                     worldData.attrib.vertices[3 * idx.vertex_index + 1],
                     worldData.attrib.vertices[3 * idx.vertex_index + 2]};
            worldData.shapesCenter[i] += p;
            worldData.shapesAABOX[i].expandBy(p);
            aabb.expandBy(p);
        }
        worldData.shapesCenter[i] /= float(shape.mesh.indices.size());
    }

    // Build BVH
    bvh = std::unique_ptr<TinyRender::AcceleratorBVH>(new TinyRender::AcceleratorBVH(this->worldData));

    const clock_t beginBVH = clock();
    bvh->build();
    std::cout << "BVH built in " << float(clock() - beginBVH) / CLOCKS_PER_SEC << "s" << std::endl;

    return true;
}

float Scene::getShapeArea(const size_t shapeID, Distribution1D& faceAreaDistribution) {
    const tinyobj::shape_t& s = worldData.shapes[shapeID];

    for (size_t i = 0; i < s.mesh.indices.size(); i += 3) {
        const int i0 = s.mesh.indices[i + 0].vertex_index;
        const int i1 = s.mesh.indices[i + 1].vertex_index;
        const int i2 = s.mesh.indices[i + 2].vertex_index;
        const v3f v0{worldData.attrib.vertices[3 * i0 + 0], worldData.attrib.vertices[3 * i0 + 1],
                     worldData.attrib.vertices[3 * i0 + 2]};
        const v3f v1{worldData.attrib.vertices[3 * i1 + 0], worldData.attrib.vertices[3 * i1 + 1],
                     worldData.attrib.vertices[3 * i1 + 2]};
        const v3f v2{worldData.attrib.vertices[3 * i2 + 0], worldData.attrib.vertices[3 * i2 + 1],
                     worldData.attrib.vertices[3 * i2 + 2]};

        const v3f e1{v1 - v0};
        const v3f e2{v2 - v0};
        const v3f e3{glm::cross(e1, e2)};
        faceAreaDistribution.add(0.5f * std::sqrt(e3.x * e3.x + e3.y * e3.y + e3.z * e3.z));
    }
    const float area = faceAreaDistribution.cdf.back();
    faceAreaDistribution.normalize();
    return area;
}

v3f Scene::getFirstLightPosition() const {
    return worldData.shapesCenter[emitters[0].shapeID];
}

v3f Scene::getFirstLightIntensity() const {
    return emitters[0].getRadiance(); // point lights are defined by intensity not radiance
}

float Scene::getShapeRadius(const size_t shapeID) const {
    assert(shapeID < worldData.shapes.size());
    v3f emitterCenter = worldData.shapesCenter[shapeID];
    return worldData.shapesAABOX[shapeID].max.x - emitterCenter.x;
}

v3f Scene::getShapeCenter(const size_t shapeID) const {
    assert(shapeID < worldData.shapes.size());
    return worldData.shapesCenter[shapeID];
}

size_t Scene::getFirstLight() const {
    if (emitters.size() <= 0) return -1;
    return emitters[0].shapeID;
}

v3f Scene::getObjectVertexPosition(size_t objectIdx, size_t vertexIdx) const {
    const tinyobj::attrib_t& sa = worldData.attrib;
    const tinyobj::shape_t& s = worldData.shapes[objectIdx];

    int idx = s.mesh.indices[vertexIdx].vertex_index;
    float x = sa.vertices[3 * idx + 0];
    float y = sa.vertices[3 * idx + 1];
    float z = sa.vertices[3 * idx + 2];
    return v3f(x,y,z);
}

v3f Scene::getObjectVertexNormal(size_t objectIdx, size_t vertexIdx) const {
    const tinyobj::attrib_t& sa = worldData.attrib;
    const tinyobj::shape_t& s = worldData.shapes[objectIdx];

    int idx_n = s.mesh.indices[vertexIdx].normal_index;
    float nx = sa.normals[3 * idx_n + 0];
    float ny = sa.normals[3 * idx_n + 1];
    float nz = sa.normals[3 * idx_n + 2];
    return glm::normalize(v3f(nx,ny,nz));
}

size_t Scene::getObjectNbVertices(size_t objectIdx) const {
    return worldData.shapes[objectIdx].mesh.indices.size();
}

int Scene::getPrimitiveID(size_t vertexIdx) const {
    return vertexIdx / 3;
}

int Scene::getMaterialID(size_t objectIdx, int primID) const {
    return worldData.shapes[objectIdx].mesh.material_ids[primID];
}

TR_NAMESPACE_END
