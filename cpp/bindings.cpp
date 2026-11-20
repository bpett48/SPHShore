// SPHShore Python extension: CPython C-API interface around the existing C++ solver.
// No third-party C++ binding dependency is required to build the wheels.
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <cstdint>
#include <cmath>
#include <limits>
#include <vector>
#include <string>

// These declarations must match the original solver C ABI.
struct WT_Config {
    int N;
    double TimeStep;
    double MinX, MinY, MinZ;
    double MaxX, MaxY, MaxZ;
    double PistonBaseX, PistonAmplitude, WaveFrequency;
    double WaveStartTime, WaveRampTime;
    double SoundSpeed, SmoothingLengthFactor;
    double DensityDiffusion, ArtificialViscosity, XSPHFactor;
    unsigned int RandomSeed;
    int WriteFiles;
    double SaveEveryTime;
};
struct WT_MapRow {
    double X,Y,Z,VX,VY,VZ;
    int Count;
};
struct WT_ParticleRow {
    int ID;
    double X,Y,Z,VX,VY,VZ,Density,Pressure;
};
struct WT_Info {
    std::uint64_t Iteration;
    double Time, Spacing, SmoothingLength, ParticleMass;
};
extern "C" {
    const char* WT_LastError();
    void* WT_Create(const WT_Config*);
    void WT_Destroy(void*);
    int WT_Step(void*, std::uint64_t);
    std::int64_t WT_GetMap(void*, WT_MapRow*, std::uint64_t);
    std::int64_t WT_GetParticles(void*, WT_ParticleRow*, std::uint64_t);
    int WT_GetInfo(void*, WT_Info*);
    int WT_SetPistonState(void*, double, double);
    int WT_UseSinusoidalPiston(void*);
}

static constexpr const char* CapsuleName = "sphshore.Simulation";
static void CapsuleDestructor(PyObject* obj) {
    void* handle = PyCapsule_GetPointer(obj, CapsuleName);
    if (handle) WT_Destroy(handle);
    else PyErr_Clear();
}
static void* Handle(PyObject* obj) {
    return PyCapsule_GetPointer(obj, CapsuleName);
}
static PyObject* NativeError() {
    const char* message = WT_LastError();
    PyErr_SetString(PyExc_RuntimeError, (message && *message) ? message : "SPHShore native solver error");
    return nullptr;
}
static PyObject* Create(PyObject*, PyObject* args) {
    PyObject* config;
    if (!PyArg_ParseTuple(args, "O", &config)) return nullptr;
    if (!PyDict_Check(config)) {
        PyErr_SetString(PyExc_TypeError, "Native config must be a dictionary");
        return nullptr;
    }
    // The Python Parameters class performs type/range checks before this call.
    auto val = [config](const char* key) -> PyObject* { return PyDict_GetItemString(config, key); };
    auto number = [&](const char* key, double& destination) -> bool {
        PyObject* o = val(key);
        if (!o) { PyErr_Format(PyExc_KeyError, "Missing config: %s", key); return false; }
        destination = PyFloat_AsDouble(o);
        if (PyErr_Occurred()) return false;
        if (!std::isfinite(destination)) { PyErr_Format(PyExc_ValueError, "%s must be finite", key); return false; }
        return true;
    };
    auto integer = [&](const char* key, long long& destination) -> bool {
        PyObject* o = val(key);
        if (!o) { PyErr_Format(PyExc_KeyError, "Missing config: %s", key); return false; }
        destination = PyLong_AsLongLong(o);
        return !PyErr_Occurred();
    };
    WT_Config c{};
    long long n=0, seed=0, write=0;
    if (!integer("N",n) || !integer("RandomSeed",seed) || !integer("WriteFiles",write)) return nullptr;
    if (n<=0 || n>std::numeric_limits<int>::max() || seed<0 || seed>UINT32_MAX || (write!=0 && write!=1)) {
        PyErr_SetString(PyExc_ValueError,"Invalid N, RandomSeed or WriteFiles"); return nullptr;
    }
    c.N=(int)n; c.RandomSeed=(unsigned int)seed; c.WriteFiles=(int)write;
    if (!number("TimeStep",c.TimeStep) ||
        !number("MinX",c.MinX) || !number("MinY",c.MinY) || !number("MinZ",c.MinZ) ||
        !number("MaxX",c.MaxX) || !number("MaxY",c.MaxY) || !number("MaxZ",c.MaxZ) ||
        !number("PistonBaseX",c.PistonBaseX) || !number("PistonAmplitude",c.PistonAmplitude) ||
        !number("WaveFrequency",c.WaveFrequency) || !number("WaveStartTime",c.WaveStartTime) ||
        !number("WaveRampTime",c.WaveRampTime) || !number("SoundSpeed",c.SoundSpeed) ||
        !number("SmoothingLengthFactor",c.SmoothingLengthFactor) ||
        !number("DensityDiffusion",c.DensityDiffusion) ||
        !number("ArtificialViscosity",c.ArtificialViscosity) ||
        !number("XSPHFactor",c.XSPHFactor) || !number("SaveEveryTime",c.SaveEveryTime)) return nullptr;
    void* handle = nullptr;
    Py_BEGIN_ALLOW_THREADS
    handle = WT_Create(&c);
    Py_END_ALLOW_THREADS
    if (!handle) return NativeError();
    PyObject* capsule = PyCapsule_New(handle, CapsuleName, CapsuleDestructor);
    if (!capsule) WT_Destroy(handle);
    return capsule;
}
static PyObject* Close(PyObject*, PyObject* args) {
    PyObject* obj;
    if (!PyArg_ParseTuple(args,"O", &obj)) return nullptr;
    void* handle=Handle(obj);
    if (!handle) return nullptr;
    // Remove destructor before removing handle, preventing a double free.
    if (PyCapsule_SetDestructor(obj,nullptr) != 0) return nullptr;
    if (PyCapsule_SetName(obj,"sphshore.ClosedSimulation") != 0) return nullptr;
    WT_Destroy(handle);
    Py_RETURN_NONE;
}
static PyObject* Step(PyObject*, PyObject* args) {
    PyObject* obj;
    unsigned long long iterations;
    if (!PyArg_ParseTuple(args,"OK", &obj, &iterations)) return nullptr;
    void* h=Handle(obj);
    if (!h) return nullptr;
    int code=0;
    Py_BEGIN_ALLOW_THREADS
    code=WT_Step(h, (std::uint64_t)iterations);
    Py_END_ALLOW_THREADS
    if (code) return NativeError();
    Py_RETURN_NONE;
}
static PyObject* Info(PyObject*, PyObject* args) {
    PyObject* obj;
    if (!PyArg_ParseTuple(args,"O", &obj)) return nullptr;
    void* h=Handle(obj); if (!h) return nullptr;
    WT_Info info{};
    if (WT_GetInfo(h,&info)) return NativeError();
    return Py_BuildValue("{s:K,s:d,s:d,s:d,s:d}",
        "iteration",(unsigned long long)info.Iteration,
        "time",info.Time, "particle_spacing",info.Spacing,
        "smoothing_length",info.SmoothingLength,"particle_mass",info.ParticleMass);
}
static PyObject* Map(PyObject*, PyObject* args) {
    PyObject* obj;
    if (!PyArg_ParseTuple(args,"O", &obj)) return nullptr;
    void* h=Handle(obj); if (!h) return nullptr;
    std::int64_t n=WT_GetMap(h,nullptr,0);
    if (n<0) return NativeError();
    if ((std::uint64_t)n > (std::uint64_t)PY_SSIZE_T_MAX) { PyErr_SetString(PyExc_OverflowError,"map too large"); return nullptr; }
    std::vector<WT_MapRow> rows((size_t)n);
    if (WT_GetMap(h,rows.data(),(std::uint64_t)n)<0) return NativeError();
    PyObject* result=PyList_New((Py_ssize_t)n);
    if (!result) return nullptr;
    for (Py_ssize_t i=0;i<(Py_ssize_t)n;i++) {
        const auto& r=rows[(size_t)i];
        PyObject* item=Py_BuildValue("(ddddddi)",r.X,r.Y,r.Z,r.VX,r.VY,r.VZ,r.Count);
        if (!item) { Py_DECREF(result); return nullptr; }
        PyList_SET_ITEM(result,i,item);
    }
    return result;
}
static PyObject* Particles(PyObject*, PyObject* args) {
    PyObject* obj;
    if (!PyArg_ParseTuple(args,"O", &obj)) return nullptr;
    void* h=Handle(obj); if (!h) return nullptr;
    std::int64_t n=WT_GetParticles(h,nullptr,0);
    if (n<0) return NativeError();
    if ((std::uint64_t)n > (std::uint64_t)PY_SSIZE_T_MAX) { PyErr_SetString(PyExc_OverflowError,"particle array too large"); return nullptr; }
    std::vector<WT_ParticleRow> rows((size_t)n);
    if (WT_GetParticles(h,rows.data(),(std::uint64_t)n)<0) return NativeError();
    PyObject* result=PyList_New((Py_ssize_t)n);
    if (!result) return nullptr;
    for (Py_ssize_t i=0;i<(Py_ssize_t)n;i++) {
        const auto& r=rows[(size_t)i];
        PyObject* item=Py_BuildValue("(idddddddd)",r.ID,r.X,r.Y,r.Z,r.VX,r.VY,r.VZ,r.Density,r.Pressure);
        if (!item) { Py_DECREF(result); return nullptr; }
        PyList_SET_ITEM(result,i,item);
    }
    return result;
}
static PyObject* SetPiston(PyObject*, PyObject* args) {
    PyObject* obj; double position,velocity;
    if (!PyArg_ParseTuple(args,"Odd", &obj,&position,&velocity)) return nullptr;
    void* h=Handle(obj); if (!h) return nullptr;
    if (WT_SetPistonState(h,position,velocity)) return NativeError();
    Py_RETURN_NONE;
}
static PyObject* Sinusoidal(PyObject*, PyObject* args) {
    PyObject* obj;
    if (!PyArg_ParseTuple(args,"O", &obj)) return nullptr;
    void* h=Handle(obj); if (!h) return nullptr;
    if (WT_UseSinusoidalPiston(h)) return NativeError();
    Py_RETURN_NONE;
}
static PyMethodDef Methods[]={
    {"create",Create,METH_VARARGS,"Create a native SPHShore simulation"},
    {"close",Close,METH_VARARGS,"Destroy a native simulation"},
    {"step",Step,METH_VARARGS,"Advance the solver by a batch of steps"},
    {"info",Info,METH_VARARGS,"Return solver state"},
    {"output_map",Map,METH_VARARGS,"Return the spatial-hash velocity map"},
    {"particles",Particles,METH_VARARGS,"Return fluid particle states"},
    {"set_piston",SetPiston,METH_VARARGS,"Set custom piston position/velocity"},
    {"use_sinusoidal_piston",Sinusoidal,METH_VARARGS,"Re-enable sinusoidal piston"},
    {nullptr,nullptr,0,nullptr}
};
static struct PyModuleDef Module={PyModuleDef_HEAD_INIT,"_core","SPHShore native C++ solver",-1,Methods};
PyMODINIT_FUNC PyInit__core(void) { return PyModule_Create(&Module); }
