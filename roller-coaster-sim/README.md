# CSCI 420 Computer Graphics, USC
## Homework 2: Roller Coaster

- **Name**: Kyuhong Lee
- **USC ID**: 7963194184
- **USC Email**: kyuhongl@usc.edu

### Build Instructions
```bash
cd hw1
make
./hw1 ../splines/rollerCoaster.sp
```

### Keyboard Controls
- `s` - Toggle skybox visibility
- `g` - Toggle ground plane visibility
- `b` - Toggle bird's eye view
- `t` - Start/stop recording frames (saves 450 frames at 15fps to frames2/)
- `ESC` - Exit program
- `SHIFT + Mouse` - Scale the scene

### Features Implemented

#### Required Features
**Levels 1-5 Complete**

#### Extra Credit Features

1. **Realistic Rail System**: 
   - **Double Rail Design**: Implemented parallel rails along the track
   - **Cylindrical Ties**: Added realistic cylindrical crossbars connecting the rails
   - **Rail Cross-Section**: Rails have a rectangular cross-section with proper lighting

2. **Skybox Environment**: 
   - **Cubemap Implementation**: Created a texture-mapped skybox using a cubemap, and the skybox was from an opensource
   - **Shader-Based Rendering**: Custom vertex and fragment shaders for skybox rendering
   - **Seamless Textures**: Properly aligned textures across cube faces
   - **Dynamic Viewing**: Skybox moves with the camera for immersive experience

3. **Bird's Eye View**: 
   - **Alternative Camera Mode**: Implemented a bird's eye view toggle for observing the entire roller coaster
   - **P.S.** It sometimes renders incorrectly. Apologies.

4. **Optimized Rendering**: 
   - **Vertex Buffer Objects**: Use of VBOs for rendering geometry
   - **Element Buffer Objects**: Used EBOs to minimize vertex duplication
   - **Shader Programs**: Specialized shader programs for different rendering needs

