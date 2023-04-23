/*
 * maximal visual madness - call me after dark
 * pixel emphasis function
 * dsp (low, high, mid)
 * starfield
 * last animation frame
 * sparkle
 * special tween non linear
 * spherical projection / off center
 * Pixel growth (starify)
 * layers
 * effect chains
 * color palette loader via png
 * zig zag function
 * game of life layer
 * lorenz attractor
 * tlfx / txmirror
 * fractalizer (barycentric)
 * mix with original
 *
 *
 *
 ffmpeg -framerate 25 -pattern_type glob -i 'walker2_anim_*.png' -i 'shona.ogg' -c:v libx264 -c:a copy -shortest -r 30 -pix_fmt yuv420p walker2_anim.mp4
 *
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <png.h>
#include <random>
#include <queue>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <fmt/core.h>
#include <vorbis/vorbisfile.h>


// using namespace std;
namespace mvm {
    const int NUMBEROFCOLORS = 1500;
    const uint8_t VALUES_18BPP[] = {0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84,
                                    88, 92, 96, 100, 104, 108, 112, 116, 120, 126, 130, 136, 140, 144, 148, 152, 156,
                                    160, 164, 168, 172, 176, 180, 184, 188, 192, 196, 200, 204, 208, 212, 216, 220, 224,
                                    228, 232, 236, 240, 244, 248, 252, 255};

/////////////////////////////////////////////////////////////////////////////
//	Floyd-Steinberg dither
/////////////////////////////////////////////////////////////////////////////

//	Floyd-Steinberg dither uses constants 7/16 5/16 3/16 and 1/16
//	But instead of using real arythmetic, I will use integer on by
//	applying shifting ( << 8 )
//	When use the constants, don't foget to shift back the result ( >> 8 )
#define    f7_16    112        //const int	f7	= (7 << 8) / 16;
#define    f5_16     80        //const int	f5	= (5 << 8) / 16;
#define    f3_16     48        //const int	f3	= (3 << 8) / 16;
#define    f1_16     16        //const int	f1	= (1 << 8) / 16;

    class Dsp {
    public:
        bool loadOgg(const std::string &filename) {
            // Open the file
            FILE *file = fopen(filename.c_str(), "rb");
            if (!file) {
                std::cerr << "Failed to open file " << filename << std::endl;
                return false;
            }

            // Initialize the Ogg Vorbis stream
            OggVorbis_File oggFile;
            if (ov_open(file, &oggFile, nullptr, 0) < 0) {
                std::cerr << "Failed to open Ogg file " << filename << std::endl;
                fclose(file);
                return false;
            }

            // Get information about the audio stream
            vorbis_info *info = ov_info(&oggFile, -1);
            channels = info->channels;
            sampleRate = info->rate;

            // Decode the audio data and store it in a vector
            totalSamples = ov_pcm_total(&oggFile, -1);

            fmt::println("DSP loading {}: {} samples/s, {} channels, {} total samples", filename, sampleRate, channels, totalSamples);
            pcmData.resize(totalSamples * channels);
            float** pcm;
            size_t startOffset = 0;
            long samplesRead = 0;
            int bitstream = 0;
            while (samplesRead < totalSamples) {
                long samplesToRead = std::min((long) totalSamples - samplesRead, 4096L);
                long bytesRead = ov_read_float(&oggFile, &pcm, samplesToRead, &bitstream);
                if (bytesRead <= 0) {
                    break;
                }

                if(bytesRead > 0)
                {
                    for (int c = 0; c < channels; ++c)
                    {
                        for (int s = 0; s < bytesRead; ++s)
                        {
                            pcmData[startOffset + s * channels + c] = pcm[c][s];
                        }
                    }
                }
                startOffset += bytesRead * channels;
                samplesRead += bytesRead;
                // pcm += bytesRead * channels;
            }
            // delete[] pcm; // TODO: check why this crashes

            // Clean up and return success
            ov_clear(&oggFile);
            fclose(file);
            return true;
        }

        void generateDefaultBands(float lowFc, float midFl, float midFh, float highFc) {
            bands.resize(3);
            bands[0].resize(totalSamples);
            bands[1].resize(totalSamples);
            bands[2].resize(totalSamples);
            std::vector<float> temp (totalSamples);
            float left, right;
            // generate low pass mono data
            fmt::println("DSP generating low pass ...");

            for (size_t s=0; s<totalSamples; s++) {
                left=pcmData[s];
                // fmt::println("fl{}", left);
                right=pcmData[s+totalSamples];
                temp[s]= (left + right) / 2.0f;
            }
            lowpassFilter(temp, bands[0], lowFc);
            // generate mid pass mono data
            fmt::println("DSP generating mid pass ...");

            for (size_t s=0; s<totalSamples; s++) {
                left=pcmData[s];
                right=pcmData[s+totalSamples];
                temp[s]= (left + right) / 2.0f;
            }
            midpassFilter(temp, bands[1], midFl, midFh);
            // generate high pass mono data
            fmt::println("DSP generating high pass ...");
            for (size_t s=0; s<totalSamples; s++) {
                left=pcmData[s];
                right=pcmData[s+totalSamples];
                temp[s]= (left + right) / 2.0f;
            }
            highpassFilter(temp, bands[2], highFc);
        }

        void generateVuMeters(uint32_t framesPerSecond) {
            size_t totalMeterSamples = std::ceil( ((double)totalSamples / sampleRate) * framesPerSecond);
            vumeters.resize(bands.size());
            for (size_t band=0; band < bands.size(); band++) {
                vumeters[band].resize(totalMeterSamples);
            }
            // we have junkSize input signals per frame
            size_t junkSize = std::ceil((double)sampleRate / framesPerSecond);
            std::vector<float> temp (junkSize);

            // generate for each frequency band
            for (size_t band=0; band < bands.size(); band++) {
                fmt::println("DSP generating VU meter data for band {} ...", band);
                size_t meterindex = 0;
                for (size_t s = 0; s < totalSamples; s += junkSize) {
                    for (size_t j = 0; j < junkSize; j++) {
                        size_t lookupIndex = s + j;
                        if (lookupIndex > totalSamples) {
                            fmt::println("Lookupindex out of bounds: {} > {}", lookupIndex, totalSamples);
                            lookupIndex = totalSamples;
                        }
                        temp[j] = bands[band][lookupIndex];
                        // fmt::println("{}",bands[band][lookupIndex] );
                    }
                    vumeters[band][meterindex] = uvMeter(temp);
                    meterindex++;
                }
            }
        }

        // Low-pass filter with cutoff frequency fc
        void lowpassFilter(const std::vector<float> &input, std::vector<float> &output, float fc) {
            float c = 1.0f / (1.0f + 2.0f * 3.14159265f * fc);
            output.resize(input.size());
            output[0] = input[0];
            for (int i = 1; i < input.size(); i++) {
                output[i] = c * (input[i] + input[i - 1]) + (1.0f - 2.0f * c) * output[i - 1];
            }
        }

        // Mid-pass filter with cutoff frequencies fl and fh
        void midpassFilter(const std::vector<float> &input, std::vector<float> &output, float fl, float fh) {
            float c1 = 2.0f * 3.14159265f * fl;
            float c2 = 2.0f * 3.14159265f * fh;
            float c = 1.0f / (1.0f + c1 / c2 + c1 * c1 / (c2 * c2));
            output.resize(input.size());
            output[0] = input[0];
            for (int i = 1; i < input.size(); i++) {
                output[i] = c * (input[i] - input[i - 1] + c1 / c2 * output[i - 1]) +
                            (1.0f - c * (1.0f + c1 / c2)) * output[i - 1];
            }
        }

        // High-pass filter with cutoff frequency fc
        void highpassFilter(const std::vector<float> &input, std::vector<float> &output, float fc) {
            float c = 1.0f / (1.0f + 2.0f * 3.14159265f * fc);
            output.resize(input.size());
            output[0] = input[0];
            for (int i = 1; i < input.size(); i++) {
                output[i] = c * (input[i] - input[i - 1]) + (1.0f - c) * output[i - 1];
            }
        }


        // Band-pass filter with center frequency fc and bandwidth bw
        void bandpassFilter(const std::vector<float> &input, std::vector<float> &output, float fc, float bw) {
            float c = 2.0f * 3.14159265f * bw;
            float a = 2.0f * c * fc;
            float b = c * c;
            float r = 1.0f / (1.0f + a + b);
            float c1 = (1.0f + a) * r;
            float c2 = -2.0f * r;
            float c3 = (1.0f - a + b) * r;
            output.resize(input.size());
            output[0] = input[0];
            output[1] = c1 * input[1] + c2 * input[0];
            for (int i = 2; i < input.size(); i++) {
                output[i] = c1 * input[i] + c2 * input[i - 1] + c3 * input[i - 2];
            }
        }

        float vuMeter(std::vector<float> &samples) {
            float sum = 0.0f;

            // Calculate the sum of the squares of each sample
            for (float sample : samples) {
                sum += sample * sample;
            }

            // Take the square root of the average of the squared samples
            float rms = std::sqrt(sum / samples.size());

            // Convert the RMS value to a dBFS value using the formula: dBFS = 20 * log10(RMS / 1.0)
            float dbfs = 20.0f * std::log10(rms / 1.0f);

            // Calculate the UV meter value by mapping the dBFS value to the appropriate range
            float vumeter = 0.0f;
            if (dbfs >= -20.0f) {
                vumeter = 0.0f;
            } else if (dbfs >= -26.0f) {
                vumeter = (dbfs + 26.0f) * 0.025f;
            } else if (dbfs >= -34.0f) {
                vumeter = (dbfs + 34.0f) * 0.025f + 0.5f;
            } else if (dbfs >= -42.0f) {
                vumeter = (dbfs + 42.0f) * 0.025f + 1.0f;
            } else if (dbfs >= -50.0f) {
                vumeter = (dbfs + 50.0f) * 0.025f + 1.5f;
            } else if (dbfs >= -58.0f) {
                vumeter = (dbfs + 58.0f) * 0.025f + 2.0f;
            } else if (dbfs >= -66.0f) {
                vumeter = (dbfs + 66.0f) * 0.025f + 2.5f;
            } else if (dbfs >= -74.0f) {
                vumeter = (dbfs + 74.0f) * 0.025f + 3.0f;
            } else if (dbfs >= -82.0f) {
                vumeter = (dbfs + 82.0f) * 0.025f + 3.5f;
            } else if (dbfs >= -90.0f) {
                vumeter = (dbfs + 90.0f) * 0.025f + 4.0f;
            } else {
                vumeter = 5.0f;
            }

            return vumeter;
        }

        float uvMeter(std::vector<float>& samples) {
            float running_sum = 0.0f;
            float running_average = 0.0f;
            float dbfs = 0.0f;
            float uvmeter = 0.0f;
            int count = 0;

            // Iterate through each sample in the vector
            for (float sample : samples) {
                // Add the sample to the running sum
                running_sum += sample;

                // Calculate the running average
                count++;
                running_average = running_sum / count;

                // Calculate the RMS value from the running average
                float rms = std::sqrt(running_average * running_average);

                // Convert the RMS value to a dBFS value
                dbfs = 20.0f * std::log10(rms / 1.0f);

                // Map the dBFS value to the appropriate range to obtain the UV meter value
                if (dbfs >= -20.0f) {
                    uvmeter = 0.0f;
                } else if (dbfs >= -26.0f) {
                    uvmeter = (dbfs + 26.0f) * 0.025f;
                } else if (dbfs >= -34.0f) {
                    uvmeter = (dbfs + 34.0f) * 0.025f + 0.5f;
                } else if (dbfs >= -42.0f) {
                    uvmeter = (dbfs + 42.0f) * 0.025f + 1.0f;
                } else if (dbfs >= -50.0f) {
                    uvmeter = (dbfs + 50.0f) * 0.025f + 1.5f;
                } else if (dbfs >= -58.0f) {
                    uvmeter = (dbfs + 58.0f) * 0.025f + 2.0f;
                } else if (dbfs >= -66.0f) {
                    uvmeter = (dbfs + 66.0f) * 0.025f + 2.5f;
                } else if (dbfs >= -74.0f) {
                    uvmeter = (dbfs + 74.0f) * 0.025f + 3.0f;
                } else if (dbfs >= -82.0f) {
                    uvmeter = (dbfs + 82.0f) * 0.025f + 3.5f;
                } else if (dbfs >= -90.0f) {
                    uvmeter = (dbfs + 90.0f) * 0.025f + 4.0f;
                } else {
                    uvmeter = 5.0f;
                }
            }

            return uvmeter;
        }




        std::vector<float> pcmData;
        int channels;
        int sampleRate;
        size_t totalSamples;

        std::vector<std::vector<float>> bands;
        std::vector<std::vector<float>> vumeters;
    };


    bool fileExists(const std::string &filename) {
        std::ifstream file(filename.c_str());
        return file.good();
    }

    struct Color {
        unsigned char r;
        unsigned char g;
        unsigned char b;
        unsigned char alpha;

        bool operator==(const Color &other) const {
            return ((r == other.r) && (g == other.g) && (b == other.b) && (alpha == other.alpha));
        }

        bool operator!=(const Color &other) const {
            return !(*this == other);
        }
    };

    inline Color mixColors(Color color1, Color color2, float weight) {
        Color result;
        result.r = (unsigned char) (color1.r * (1 - weight) + color2.r * weight);
        result.g = (unsigned char) (color1.g * (1 - weight) + color2.g * weight);
        result.b = (unsigned char) (color1.b * (1 - weight) + color2.b * weight);
        result.alpha = (unsigned char) (color1.alpha * (1 - weight) + color2.alpha * weight);
        return result;
    }

    inline Color addColors(Color color1, Color color2) {
        Color result;
        result.r = std::clamp(color1.r + color2.r, 0, 255);
        result.g = std::clamp(color1.g + color2.g, 0, 255);
        result.b = std::clamp(color1.b + color2.b, 0, 255);
        result.alpha = std::clamp(color1.alpha + color2.alpha, 0, 255);
        // fmt::print("{} {} {}\n", color1.r, color2.r, result.r);
        return result;
    }

// Define function to interpolate between two RGBA colors
    Color interpolate(Color c1, Color c2, float fraction) {
        float r = c1.r * (1 - fraction) + c2.r * fraction;
        float g = c1.g * (1 - fraction) + c2.g * fraction;
        float b = c1.b * (1 - fraction) + c2.b * fraction;
        float alpha = c1.alpha * (1 - fraction) + c2.alpha * fraction;

        Color result = {result.r = round(r), result.g = round(g), result.b = round(b), result.alpha = round(alpha)};
        return result;
    }


    struct Pixel {
        uint32_t x;
        uint32_t y;
        Color color;
    };

// Define function to interpolate between two pixels
    Pixel interpolate(const Pixel &p1, const Pixel &p2, float fraction) {
        float x = p1.x * (1 - fraction) + p2.x * fraction;
        float y = p1.y * (1 - fraction) + p2.y * fraction;
        Pixel result;
        result.x = round(x);
        result.y = round(y);
        result.color = interpolate(p1.color, p2.color, fraction);
        return result;
    }

    std::vector<Color> generate_heatmap_color_table(int n) {
        std::vector<Color> color_table;
        for (int i = 0; i < n; i++) {
            double t = (double) i / (double) (n - 1);
            unsigned char r = round(255 * pow(1.0 - t, 2.0));
            unsigned char g = round(255 * 2.0 * (1.0 - t) * t);
            unsigned char b = round(255 * pow(t, 2.0));
            Color color = {r, g, b, 255};
            color_table.push_back(color);
        }
        return color_table;
    }

    class Image {
    public:
        Image(int width, int height) {
            this->width = width;
            this->height = height;
            data.resize(width * height * 4);
            this->heatmapColors = generate_heatmap_color_table(NUMBEROFCOLORS);
        }

        Image(int width, int height, Color colorBackground) {
            this->width = width;
            this->height = height;
            data.resize(width * height * 4);
            this->heatmapColors = generate_heatmap_color_table(NUMBEROFCOLORS);
            this->fill(colorBackground);
        }

        Image(std::string filename) {
            width = 0;
            height = 0;
            FILE *fp = fopen(filename.c_str(), "rb");
            if (!fp) {
                throw std::runtime_error("Error opening PNG file.");
            }

            png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
            if (!png_ptr) {
                fclose(fp);
                throw std::runtime_error("Error creating PNG read structure.");
            }

            png_infop info_ptr = png_create_info_struct(png_ptr);
            if (!info_ptr) {
                png_destroy_read_struct(&png_ptr, NULL, NULL);
                fclose(fp);
                throw std::runtime_error("Error creating PNG info structure.");
            }

            if (setjmp(png_jmpbuf(png_ptr))) {
                png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
                fclose(fp);
                throw std::runtime_error("Error reading PNG file.");
            }

            png_init_io(png_ptr, fp);
            png_read_info(png_ptr, info_ptr);

            width = png_get_image_width(png_ptr, info_ptr);
            height = png_get_image_height(png_ptr, info_ptr);
            png_byte color_type = png_get_color_type(png_ptr, info_ptr);
            png_byte bit_depth = png_get_bit_depth(png_ptr, info_ptr);

            if (color_type != PNG_COLOR_TYPE_RGBA || bit_depth != 8) {
                png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
                fclose(fp);
                throw std::runtime_error("Error loading PNG: unsupported color type or bit depth.");
            }
            png_bytep *row_pointers = (png_bytep *) malloc(sizeof(png_bytep) * height);
            for (png_uint_32 y = 0; y < height; y++) {
                row_pointers[y] = (png_byte *) malloc(png_get_rowbytes(png_ptr, info_ptr));
            }

            png_read_image(png_ptr, row_pointers);
            data.reserve(width * height * 4);
            for (png_uint_32 y = 0; y < height; y++) {
                png_bytep row = row_pointers[y];
                for (png_uint_32 x = 0; x < width; x++) {
                    png_bytep px = &(row[x * 4]);
                    data.push_back(px[0]);
                    data.push_back(px[1]);
                    data.push_back(px[2]);
                    data.push_back(px[3]);
                }
            }

            for (png_uint_32 y = 0; y < height; y++) {
                free(row_pointers[y]);
            }
            free(row_pointers);

            png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
            fclose(fp);
        }

        bool savePNG(std::string filename) {
            FILE *fp = fopen(filename.c_str(), "wb");
            if (!fp) {
                return false;
            }

            png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
            if (!png_ptr) {
                fclose(fp);
                return false;
            }

            png_infop info_ptr = png_create_info_struct(png_ptr);
            if (!info_ptr) {
                png_destroy_write_struct(&png_ptr, NULL);
                fclose(fp);
                return false;
            }

            if (setjmp(png_jmpbuf(png_ptr))) {
                png_destroy_write_struct(&png_ptr, &info_ptr);
                fclose(fp);
                return false;
            }

            png_init_io(png_ptr, fp);

            png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
                         PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

            png_write_info(png_ptr, info_ptr);

            std::vector<png_bytep> row_pointers(height);
            for (int y = 0; y < height; ++y) {
                row_pointers[y] = (png_bytep) &data[y * width * 4];
            }

            png_write_image(png_ptr, row_pointers.data());
            png_write_end(png_ptr, NULL);

            png_destroy_write_struct(&png_ptr, &info_ptr);
            fclose(fp);
            return true;
        }

        void setAlpha(unsigned char alpha) {
            for (int i = 3; i < this->data.size(); i += 4) {
                this->data.at(i) = alpha;
            }
        }

        void fill(Color color) {
            for (int i = 0; i < this->data.size(); i += 4) {
                this->data.at(i) = color.r;
                this->data.at(i + 1) = color.g;
                this->data.at(i + 2) = color.b;
                this->data.at(i + 3) = color.alpha;
            }
        }

        // Get the index of a pixel in the image vector given its (x, y) coordinates
        [[nodiscard]] inline size_t pixelIndex(size_t x, size_t y) const {
            size_t index = (y * width + x) * 4;
            if (index > data.size())
                fmt::println("pixelIndex out range at {},{}: index {} > {}", x, y, index, data.size());
            return index;
        }

        uint32_t width = 0;
        uint32_t height = 0;
        std::vector<unsigned char> data;
        std::vector<Pixel> undoBuffer;
        std::vector<Color> heatmapColors;
    };

    inline bool getPixel(Image &image, int x, int y, Color &color) {
        if (x < 0 || x >= image.width || y < 0 || y >= image.height) return false; // Out of bounds
        int index = 4 * (y * image.width + x);
        color.r = image.data[index];
        color.g = image.data[index + 1];
        color.b = image.data[index + 2];
        color.alpha = image.data[index + 3];
        return true;
    }

    inline Color getPixel(Image &image, int x, int y) {
        int index = 4 * (y * image.width + x);
        Color color = {image.data[index], image.data[index + 1], image.data[index + 2], image.data[index + 3]};
        return color;
    }

    inline Color getPixel(Image &image, uint32_t index) {
        Color color;
        color.r = image.data[index];
        color.g = image.data[index + 1];
        color.b = image.data[index + 2];
        color.alpha = image.data[index + 3];
        return color;
    }

// Function to plot a pixel at a given x and y position with a given color
    inline void setPixel(Image &image, int x, int y, Color color) {
        int index = 4 * (y * image.width + x);
        image.data[index] = color.r;
        image.data[index + 1] = color.g;
        image.data[index + 2] = color.b;
        image.data[index + 3] = color.alpha;
    }

    inline void setPixel(Image &image, int x, int y, Color colorNew, float weight) {
        int index = 4 * (y * image.width + x);
        Color colorCurrent = getPixel(image, index);
        Color colorMixed = mixColors(colorCurrent, colorNew, weight);
        image.data[index] = colorMixed.r;
        image.data[index + 1] = colorMixed.g;
        image.data[index + 2] = colorMixed.b;
        image.data[index + 3] = colorMixed.alpha;
    }

    inline void setPixelAdd(Image &image, int x, int y, Color colorNew) {
        int index = 4 * (y * image.width + x);
        Color colorCurrent = getPixel(image, index);
        Color colorMixed = addColors(colorCurrent, colorNew);
        image.data[index] = colorMixed.r;
        image.data[index + 1] = colorMixed.g;
        image.data[index + 2] = colorMixed.b;
        image.data[index + 3] = colorMixed.alpha;
    }


    void gaussianBlurRGBA(Image &image) {
        // Create a 3x3 Gaussian kernel
        std::vector<float> kernel = {1.f / 16.f, 2.f / 16.f, 1.f / 16.f,
                                     2.f / 16.f, 4.f / 16.f, 2.f / 16.f,
                                     1.f / 16.f, 2.f / 16.f, 1.f / 16.f};

        // Iterate over each pixel in the image
        for (uint32_t y = 1; y < image.height - 1; ++y) {
            for (uint32_t x = 1; x < image.width - 1; ++x) {
                // Calculate the weighted sum of neighboring pixels
                float sumR = 0.f, sumG = 0.f, sumB = 0.f, sumA = 0.f;
                for (int j = -1; j <= 1; ++j) {
                    for (int i = -1; i <= 1; ++i) {
                        uint32_t index = ((y + j) * image.width + x + i) * 4;
                        float weight = kernel[(j + 1) * 3 + i + 1];
                        sumR += weight * static_cast<float>(image.data[index]);
                        sumG += weight * static_cast<float>(image.data[index + 1]);
                        sumB += weight * static_cast<float>(image.data[index + 2]);
                        sumA += weight * static_cast<float>(image.data[index + 3]);
                    }
                }

                // Set the pixel to the average of the neighboring pixels
                uint32_t index = (y * image.width + x) * 4;
                image.data[index] = static_cast<unsigned char>(std::round(sumR));
                image.data[index + 1] = static_cast<unsigned char>(std::round(sumG));
                image.data[index + 2] = static_cast<unsigned char>(std::round(sumB));
                image.data[index + 3] = static_cast<unsigned char>(std::round(sumA));
            }
        }
    }

    void replaceColor(Image &image, Color colorFind, Color colorReplace) {
        // Iterate over each pixel in the image
        for (uint32_t y = 1; y < image.height - 1; ++y) {
            for (uint32_t x = 1; x < image.width - 1; ++x) {
                if (getPixel(image, x, y) == colorFind) {
                    setPixel(image, x, y, colorReplace);
                }
            }
        }
    }


//	Color Floyd-Steinberg dither using 18 bits per pixel (6 bit per color plane)
//  https://www.codeproject.com/Articles/5259216/Dither-Ordered-and-Floyd-Steinberg-Monochrome-Colo
    void makeDitherFSRgb18bpp(Image &image) noexcept {
        const int size = image.width * image.height;

        int *errorB = (int *) malloc(size * sizeof(int));
        int *errorG = (int *) malloc(size * sizeof(int));
        int *errorR = (int *) malloc(size * sizeof(int));

        //	Clear the errors buffer.
        memset(errorB, 0, size * sizeof(int));
        memset(errorG, 0, size * sizeof(int));
        memset(errorR, 0, size * sizeof(int));

        //~~~~~~~~

        int i = 0;

        for (int y = 0; y < image.height; y++) {
            uint8_t *prow = &image.data[y * image.width * 4];

            for (int x = 0; x < image.width; x++, i++) {
                const int blue = prow[x * 4 + 0];
                const int green = prow[x * 4 + 1];
                const int red = prow[x * 4 + 2];

                int newValB = (int) blue + (errorB[i] >> 8);    //	PixelBlue  + error correctionB
                int newValG = (int) green + (errorG[i] >> 8);    //	PixelGreen + error correctionG
                int newValR = (int) red + (errorR[i] >> 8);    //	PixelRed   + error correctionR

                //	The error could produce values beyond the borders, so need to keep the color in range
                int idxR = std::clamp(newValR, 0, 255);
                int idxG = std::clamp(newValG, 0, 255);
                int idxB = std::clamp(newValB, 0, 255);

                int newcR = VALUES_18BPP[idxR >> 2];    //	x >> 2 is the same as x / 4
                int newcG = VALUES_18BPP[idxG >> 2];    //	x >> 2 is the same as x / 4
                int newcB = VALUES_18BPP[idxB >> 2];    //	x >> 2 is the same as x / 4

                prow[x * 4 + 0] = newcB;
                prow[x * 4 + 1] = newcG;
                prow[x * 4 + 2] = newcR;

                int cerrorB = newValB - newcB;
                int cerrorG = newValG - newcG;
                int cerrorR = newValR - newcR;

                int idx = i + 1;
                if (x + 1 < image.width) {
                    errorR[idx] += (cerrorR * f7_16);
                    errorG[idx] += (cerrorG * f7_16);
                    errorB[idx] += (cerrorB * f7_16);
                }

                idx += image.width - 2;
                if (x - 1 > 0 && y + 1 < image.height) {
                    errorR[idx] += (cerrorR * f3_16);
                    errorG[idx] += (cerrorG * f3_16);
                    errorB[idx] += (cerrorB * f3_16);
                }

                idx++;
                if (y + 1 < image.height) {
                    errorR[idx] += (cerrorR * f5_16);
                    errorG[idx] += (cerrorG * f5_16);
                    errorB[idx] += (cerrorB * f5_16);
                }

                idx++;
                if (x + 1 < image.width && y + 1 < image.height) {
                    errorR[idx] += (cerrorR * f1_16);
                    errorG[idx] += (cerrorG * f1_16);
                    errorB[idx] += (cerrorB * f1_16);
                }
            }
        }

        free(errorB);
        free(errorG);
        free(errorR);
    }

    bool hasNeighborOfColorC(Image &image, int x, int y, Color color) {
        // Check each neighboring pixel
        Color colorCurrent;
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                int neighborX = x + i;
                int neighborY = y + j;

                // Check if neighbor is within image bounds
                if (neighborX >= 0 && neighborX < image.width && neighborY >= 0 && neighborY < image.height) {
                    // Check if neighbor has colorC
                    if (getPixel(image, neighborX, neighborY, colorCurrent)) {
                        if (colorCurrent == color) return true;
                    }
                }
            }
        }

        // No neighboring pixel has colorC
        return false;
    }

    size_t edgeFill(Image &image, int x, int y, Color newColor, Color oldColor, Color edgeColor) {
        size_t pixelsFilled = 0;
        // std::cout << "ff " << x << "," << y << " " << pixelsFilled << "\n";
        // Base cases
        Color color;
        std::queue<std::pair<int, int>> pixels;
        pixels.push({x, y});

        while (!pixels.empty()) {
            std::pair<int, int> currPixel = pixels.front();
            pixels.pop();
            int i = currPixel.first;
            int j = currPixel.second;

            if (!getPixel(image, i, j, color)) {
                std::cout << "filled\n";
                continue;
            }

            if (color != oldColor) { continue; }

            pixelsFilled++;
            // newColor = image.heatmapColors[pixelsFilled % NUMBEROFCOLORS];
            setPixel(image, i, j, newColor);
            // std::cout << "ff2 " << i << "," << j << " " << pixelsFilled << "\n";
            if (!hasNeighborOfColorC(image, i, j, edgeColor)) {
                Pixel pixel = {(uint32_t) i, (uint32_t) j, color};
                image.undoBuffer.push_back(pixel);
            }

            // Check neighboring pixels for the old color before pushing them onto the stack
            if ((getPixel(image, i + 1, j, color)) && (color == oldColor)) pixels.push({i + 1, j});
            if ((getPixel(image, i - 1, j, color)) && (color == oldColor)) pixels.push({i - 1, j});

            if ((getPixel(image, i, j + 1, color)) && (color == oldColor)) pixels.push({i, j + 1});
            if ((getPixel(image, i, j - 1, color)) && (color == oldColor)) pixels.push({i, j - 1});

        }

        // std::cout << "ffe " << x << "," << y << " " << pixelsFilled << "\n";
        return pixelsFilled;
    }

    size_t edgeFill(Image &image, int x, int y, Color newColor, Color edgeColor, bool verbose) {
        image.undoBuffer.clear();
        image.undoBuffer.reserve(image.width * image.height);
        if (verbose) std::cout << "Edge Fill: " << x << "," << y << "\n";
        Color oldColor;
        size_t pixelsFilled = 0;
        if (getPixel(image, x, y, oldColor)) {
            if ((oldColor != newColor) && (oldColor != edgeColor)) {
                pixelsFilled = edgeFill(image, x, y, newColor, oldColor, edgeColor);
                if (verbose) {
                    auto edgePixels = pixelsFilled - image.undoBuffer.size();
                    std::cout << "Reverting non neighbours at (" << x << "," << y << ") :" << image.undoBuffer.size()
                              << "\n";
                    std::cout << "Edge pixels:" << edgePixels << "\n";
                }
                for (const auto &pixel: image.undoBuffer) {
                    setPixel(image, pixel.x, pixel.y, pixel.color);
                }

            } else {
                if (verbose) std::cout << "Fill skipped at (" << x << "," << y << ")\n";
            }
        }

        return pixelsFilled;
    }


    size_t
    floodFill(Image &image, int x, int y, Color newColor, Color oldColor, const bool test, const size_t maxPixels) {
        size_t pixelsFilled = 0;
        // std::cout << "ff " << x << "," << y << " " << pixelsFilled << "\n";
        // Base cases
        Color color;
        std::queue<std::pair<int, int>> pixels;
        pixels.push({x, y});

        while (!pixels.empty()) {
            std::pair<int, int> currPixel = pixels.front();
            pixels.pop();
            int i = currPixel.first;
            int j = currPixel.second;

            if (!getPixel(image, i, j, color)) {
                std::cout << "filled\n";
                continue;
            }

            if (color != oldColor) { continue; }

            pixelsFilled++;
            // newColor = image.heatmapColors[pixelsFilled % NUMBEROFCOLORS];
            setPixel(image, i, j, newColor);
            // std::cout << "ff2 " << i << "," << j << " " << pixelsFilled << "\n";
            if (test) {
                Pixel pixel = {(uint32_t) i, (uint32_t) j, color};
                image.undoBuffer.push_back(pixel);
                if (pixelsFilled >= maxPixels) return pixelsFilled;
            }

            // Check neighboring pixels for the old color before pushing them onto the stack
            if ((getPixel(image, i + 1, j, color)) && (color == oldColor)) pixels.push({i + 1, j});
            if ((getPixel(image, i - 1, j, color)) && (color == oldColor)) pixels.push({i - 1, j});

            if ((getPixel(image, i, j + 1, color)) && (color == oldColor)) pixels.push({i, j + 1});
            if ((getPixel(image, i, j - 1, color)) && (color == oldColor)) pixels.push({i, j - 1});

        }

        // std::cout << "ffe " << x << "," << y << " " << pixelsFilled << "\n";
        return pixelsFilled;
    }

    size_t floodFill(Image &image, int x, int y, Color newColor, Color excludeColor, bool test, size_t maxPixels,
                     bool verbose) {
        image.undoBuffer.clear();
        image.undoBuffer.reserve(maxPixels);
        if (verbose) std::cout << "Flood Fill: " << x << "," << y << "\n";
        Color oldColor;
        size_t pixelsFilled = 0;
        if (getPixel(image, x, y, oldColor)) {
            if ((oldColor != newColor) && (oldColor != excludeColor)) {
                pixelsFilled = floodFill(image, x, y, newColor, oldColor, test, maxPixels);
                if (test && (pixelsFilled >= maxPixels)) {
                    // revert fill operation
                    if (verbose)
                        std::cout << "Reverting fill at (" << x << "," << y << ") :" << maxPixels << " / "
                                  << image.undoBuffer.size() << "\n";
                    for (const auto &pixel: image.undoBuffer) {
                        setPixel(image, pixel.x, pixel.y, pixel.color);
                    }
                } else {
                    if (verbose) std::cout << "Filled at (" << x << "," << y << ") :" << pixelsFilled << "\n";
                }

            } else {
                if (verbose) std::cout << "Fill skipped at (" << x << "," << y << ")\n";
            }
        }

        return pixelsFilled;
    }

    void drawLine(Image &image, int x1, int y1, int x2, int y2, Color color, int threshold_manhattan_length = 0) {

        int dx = abs(x2 - x1);
        int dy = abs(y2 - y1);
        int sx = (x1 < x2) ? 1 : -1;
        int sy = (y1 < y2) ? 1 : -1;
        int err = dx - dy;

        int manhattanDistance = (dx + dy) / 2;

        if (threshold_manhattan_length) {
            if (manhattanDistance > threshold_manhattan_length) return;
        }


        while (true) {
            setPixel(image, x1, y1, color);

            if (x1 == x2 && y1 == y2) {
                break;
            }

            int e2 = 2 * err;

            if (e2 > -dy) {
                err -= dy;
                x1 += sx;
            }

            if (e2 < dx) {
                err += dx;
                y1 += sy;
            }
        }
    }

    void drawLinesAroundCenter(Image &image, int x1, int y1, int x2, int y2, Color color,
                               int threshold_manhattan_length = 0) {

        int center_x = image.width / 2;
        int center_y = image.height / 2;

        int x1m, x2m, y1m, y2m;
        int x1_offset = x1 - center_x;
        int y1_offset = y1 - center_y;
        int x2_offset = x2 - center_x;
        int y2_offset = y2 - center_y;
        if (x1_offset < 0) x1m = center_x + abs(x1_offset); else x1m = center_x - x1_offset;
        if (y1_offset < 0) y1m = center_y + abs(y1_offset); else y1m = center_y - y1_offset;
        if (x2_offset < 0) x2m = center_x + abs(x2_offset); else x2m = center_x - x2_offset;
        if (y2_offset < 0) y2m = center_y + abs(y2_offset); else y2m = center_y - y2_offset;


        drawLine(image, x1, y1, x2, y2, color, threshold_manhattan_length);
        drawLine(image, x1, y1m, x2, y2m, color, threshold_manhattan_length);
        drawLine(image, x1m, y1, x2m, y2, color, threshold_manhattan_length);
        drawLine(image, x1m, y1m, x2m, y2m, color, threshold_manhattan_length);
    }

    std::vector<Pixel> displacePixels(Image &image, int n) {
        std::cout << "displacing pixels by radius " << n << "\n";
        std::random_device rd;  // obtain a random seed from hardware
        std::mt19937 gen(rd());  // seed the generator
        std::uniform_int_distribution<> distr(-n, n);  // define the range
        std::vector<Pixel> pdata;
        pdata.reserve(image.width * image.height);
        Color color;
        for (int y = 0; y < image.height; y++) {
            for (int x = 0; x < image.width; x++) {
                int x_offset = distr(gen);
                int y_offset = distr(gen);
                int x_new = x + x_offset;
                int y_new = y + y_offset;
                if (x_new > 0 && x_new < image.width && y_new > 0 && y_new < image.height) {
                    getPixel(image, x, y, color);
                    Pixel pixel = {(uint32_t) x_new, (uint32_t) y_new, color};
                    pdata.emplace_back(pixel);
                }
            }
        }
        return pdata;
    }

    size_t getPixelsNotHavingBackground(Image &image, Color colorBackground, std::vector<Pixel> &pixels) {
        Color color;
        size_t count = 0;
        for (uint32_t y = 0; y < image.height; y++) {
            for (uint32_t x = 0; x < image.width; x++) {
                getPixel(image, x, y, color);
                if (color != colorBackground) {
                    Pixel pixel = {x, y, color};
                    pixels.emplace_back(pixel);
                    count++;
                }
            }
        }
        return count;
    }

    void drawPixelsIfNotProtected(Image &image, std::vector<Pixel> &pixels, std::vector<Color> colorsProtected,
                                  float weight) {
        for (const auto &pixel: pixels) {
            Color currentColor = getPixel(image, pixel.x, pixel.y);
            bool pixelProtected = false;
            for (const auto &colorCheck: colorsProtected) {
                if (currentColor == colorCheck) pixelProtected = true;
            }
            if (!pixelProtected) setPixel(image, pixel.x, pixel.y, pixel.color, weight);
        }
    }

#include <cmath>
#include <vector>


// Rotate an image represented by a vector of bytes by a given angle (in degrees)
    void rotateImage(Image &image, Color colorBackground, float degrees) {
        // Convert the angle to radians
        float radians = degrees * M_PI / 180.0f;
        float cosTheta = std::cos(radians);
        float sinTheta = std::sin(radians);

        // Get the dimensions of the image
        size_t width = image.width;
        size_t height = image.height;

        // Find the center of the image
        float centerX = static_cast<float>(width) / 2.0f;
        float centerY = static_cast<float>(height) / 2.0f;

        // Create a new image buffer to store the rotated image
        Image rotatedImage(image.width, image.height, colorBackground);
        // rotatedImage.fill(colorBackground);
        // Loop over each pixel in the output image
        for (size_t y = 0; y < height; y++) {
            for (size_t x = 0; x < width; x++) {
                // Compute the coordinates of the pixel in the input image
                float newX = cosTheta * (x - centerX) + sinTheta * (y - centerY) + centerX;
                float newY = -sinTheta * (x - centerX) + cosTheta * (y - centerY) + centerY;

                // Get the pixel value at the new coordinates using bilinear interpolation
                if (newX >= 0 && newX < width - 1 && newY >= 0 && newY < height - 1) {
                    size_t x0 = static_cast<size_t>(newX);
                    size_t y0 = static_cast<size_t>(newY);
                    size_t x1 = x0 + 1;
                    size_t y1 = y0 + 1;
                    float alpha = newX - x0;
                    float beta = newY - y0;
                    unsigned char *pixel00 = &image.data[image.pixelIndex(x0, y0)];
                    unsigned char *pixel01 = &image.data[image.pixelIndex(x0, y1)];
                    unsigned char *pixel10 = &image.data[image.pixelIndex(x1, y0)];
                    unsigned char *pixel11 = &image.data[image.pixelIndex(x1, y1)];
                    for (size_t c = 0; c < 4; c++) {
                        float top =
                                static_cast<float>(pixel00[c]) * (1 - alpha) + static_cast<float>(pixel10[c]) * alpha;
                        float bottom =
                                static_cast<float>(pixel01[c]) * (1 - alpha) + static_cast<float>(pixel11[c]) * alpha;
                        float value = top * (1 - beta) + bottom * beta;
                        auto targetPixelIndex = image.pixelIndex(x, y) + c;
                        if (targetPixelIndex < image.data.size())
                            rotatedImage.data[targetPixelIndex] = static_cast<unsigned char>(value);
                    }
                }
            }
        }

        // Copy the rotated image data back to the input vector
        image.data = std::move(rotatedImage.data);
    }


    std::vector<unsigned char> random_walker(uint64_t n) {

        std::random_device rd;  // obtain a random seed from hardware
        std::mt19937 gen(rd());  // seed the generator
        std::uniform_int_distribution<> distr(0, 7);  // define the range
        std::vector<unsigned char> pdata(n);
        for (size_t i = 0; i < n; i++) {
            pdata[i] = distr(gen);
            // std::cout << pdata[i]<< "\n";
        }
        std::cout << "random walker: " << pdata.size() << "\n";
        return pdata;
    }


    std::vector<unsigned char> fib(int numbers, double fraction) {
        std::vector<unsigned char> pdata(numbers);
        int count = 0;
        int sum = 0;
        int index = 0;

        for (int i = 0; i < numbers; i++) {
            sum += i * fraction;
            pdata[i] = sum % 255;
        }

        return pdata;
    }

    int turtle_square(int x, int y, int stepsize, Image &image, int maxiterations, Color color, int threshold,
                      double fibfraction) {
        auto numbers = fib(maxiterations, fibfraction);
        int xo, yo = 0;
        for (auto &number: numbers) {
            int direction = number % 4;
            // std::cout << (int) number << " " << direction << "\n";
            switch (direction) {
                case 0:
                    xo = 0;
                    yo = -1 * stepsize;
                    break;
                case 1:
                    xo = 1 * stepsize;
                    yo = 0;
                    break;
                case 2:
                    xo = 0;
                    yo = 1 * stepsize;
                    break;
                case 3:
                    xo = -1 * stepsize;
                    yo = 0;
                    break;
            }
            int x2 = abs(int((x + xo) % image.width));
            int y2 = abs(int((y + yo) % image.height));
            drawLine(image, x, y, x2, y2, color, threshold);
            x = x2;
            y = y2;
        }
        return 0;
    }

    int
    turtle(int x, int y, int stepsize, Image &image, int maxiterations, Color color, int threshold, double fibfraction,
           int drawingMode) {
        auto numbers = fib(maxiterations, fibfraction);

        int xo, yo = 0;
        int alpha = 0;
        size_t index = 0;
        for (auto &number: numbers) {
            int direction = number % 8;
            // std::cout << "d: " << direction << "\n";
            switch (direction) {
                case 0:
                    xo = 0;
                    yo = -stepsize;
                    break;
                case 1:
                    xo = stepsize;
                    yo = -stepsize;
                    break;
                case 2:
                    xo = stepsize;
                    yo = 0;
                    break;
                case 3:
                    xo = stepsize;
                    yo = stepsize;
                    break;
                case 4:
                    xo = 0;
                    yo = stepsize;
                    break;
                case 5:
                    xo = -stepsize;
                    yo = stepsize;
                    break;
                case 6:
                    xo = -stepsize;
                    yo = 0;
                    break;
                case 7:
                    xo = -stepsize;
                    yo = -stepsize;
                    break;
            }

            int x2 = abs(int((x + xo) % image.width));
            int y2 = abs(int((y + yo) % image.height));
            // color = colors[(index % NUMBEROFCOLORS)];
            switch (drawingMode) {
                case 0:
                    drawLine(image, x, y, x2, y2, color, threshold);
                    break;
                case 1:
                    drawLinesAroundCenter(image, x, y, x2, y2, color, threshold);
                    break;
            }
            x = x2;
            y = y2;

            index++;
        }
        return 0;
    }

    std::string walker(double fibFraction, bool regenerate = false) {
        std::cout << "walker 2.0\n";

        const int stepSize = 10;
        const int iterations = 1000;
        const int fillThreshold = 1200;
        const bool verbose = false;
        const int drawingMode = 1;

        std::string filename = "walker2_";
        filename += std::to_string(fibFraction);
        filename += "_" + std::to_string(iterations);
        filename += "_" + std::to_string(fillThreshold);
        filename += "D" + std::to_string(drawingMode);
        filename += ".png";

        if (fileExists(filename) && !(regenerate)) {
            fmt::println("{} exists - skipping creature generation", filename);
            return filename;
        }
        Image image(1000, 1000);
        std::cout << "allocated " << image.data.size() / (1024 * 1024) << "MB Ram for Image\n";
        std::cout << "generating RGBA image\n";

        Color color = {255, 255, 255, 255};
        Color excludeColor = color;
        std::cout << "turtle ....\n";
        std::cout << "width: " << image.width << "\n";
        std::cout << "height: " << image.height << "\n";
        std::cout << "fraction: " << fibFraction << "\n";
        std::cout << "step size: " << stepSize << "\n";

        turtle(image.width / 2, image.height / 2, stepSize, image, iterations, color, stepSize * 2.5, fibFraction,
               drawingMode);
        // turtle_square(image.width / 2, image.height / 2, stepSize, image, 5000, color, stepSize*2.5, fibFraction);

        std::cout << "alpha 100% for all pixels ....\n";
        image.setAlpha(255);

        color = {255, 0, 0, 255};
        std::vector<Color> colorTable(9);
        colorTable[0] = {144, 41, 27, 255};
        colorTable[1] = {213, 78, 49, 255};
        colorTable[2] = {235, 96, 67, 255};
        colorTable[3] = {39, 89, 84, 255};
        colorTable[4] = {62, 135, 129, 255};
        colorTable[5] = {90, 191, 180, 255};
        colorTable[6] = {243, 63, 51, 255};
        colorTable[7] = {178, 163, 40, 255};
        colorTable[8] = {236, 218, 66, 255};
        std::cout << "filling magic ....\n";

        std::random_device rd;  // obtain a random seed from hardware
        std::mt19937 gen(0);  // seed the generator
        std::uniform_int_distribution<> distr(0, colorTable.size() - 1);  // define the range

        int verticalFillIndex = 0;
        for (int y = 0; y < image.height; y += stepSize) {
            verticalFillIndex++;
            color = colorTable[verticalFillIndex % colorTable.size()];
            for (int x = 0; x < image.width; x += stepSize) {
                floodFill(image, x + 2, y + 2, color, excludeColor, true, 1200, verbose);
            }
        }
        std::cout << "edge filling magic ....\n";

        Color edgeColor = {255, 255, 255, 255};
        edgeFill(image, 0, 0, color, edgeColor, verbose);
        color = {255, 255, 0, 255};
        edgeColor = {255, 0, 0, 255};
        edgeFill(image, 0, 0, color, edgeColor, verbose);
        color = {0, 0, 0, 254};
        edgeColor = {255, 255, 0, 255};
        edgeFill(image, 0, 0, color, edgeColor, verbose);
        edgeColor = {0, 0, 0, 254};
        color = {255, 255, 0, 255};
        edgeFill(image, 0, 0, color, edgeColor, verbose);


        // ----- pixel displacement -------------------------------------
        std::vector<Color> colorsProtected;
        colorsProtected.push_back(colorTable[0]);
        colorsProtected.push_back(colorTable[1]);
        colorsProtected.push_back(colorTable[2]);
        colorsProtected.push_back(colorTable[3]);
        colorsProtected.push_back(colorTable[4]);
        colorsProtected.push_back(colorTable[5]);
        colorsProtected.push_back(colorTable[6]);
        colorsProtected.push_back(colorTable[7]);
        colorsProtected.push_back(colorTable[8]);
        colorsProtected.push_back({255, 255, 255, 255});
        colorsProtected.push_back({255, 255, 0, 255});

        auto pixels = displacePixels(image, 5);
        drawPixelsIfNotProtected(image, pixels, colorsProtected, 0.8f);

        pixels = displacePixels(image, 20);
        drawPixelsIfNotProtected(image, pixels, colorsProtected, 0.6f);

        pixels = displacePixels(image, 40);
        drawPixelsIfNotProtected(image, pixels, colorsProtected, 0.4f);

        // replace yellow with black
        Color colorFind = {255, 255, 0, 255};
        Color colorReplace = {0, 0, 0, 255};
        replaceColor(image, colorFind, colorReplace);


        // std::cout << "gaussian blur 3x3 ....\n";
        // gaussianBlurRGBA(image);

        color = {255, 255, 255, 255};
        turtle(image.width / 2, image.height / 2, stepSize, image, iterations, color, stepSize * 2.5, fibFraction,
               drawingMode);

        // std::cout << "dithering magic ....\n";
        // Perform dithering on the image
        // makeDitherFSRgb18bpp(image);

        std::cout << "writing file\n";


        // Write the image to a file in PNG format
        if (image.savePNG(filename)) {
            std::cout << "Das Bild " << filename << " wurde erfolgreich gespeichert." << std::endl;
        } else {
            std::cout << "Fehler beim Speichern des Bildes: " << filename << std::endl;
        }
        return filename;
    }

    int
    interpolate(Dsp &dsp, std::string file1, std::string file2, const std::string fileTarget, size_t frames, size_t startFrame,
                Color colorBackground, bool rotate, bool shuffleStartingPositions = false) {
        // walker();

        try {
            Image image1(file1);
            Image image2(file2);
            // get a list of all pixels not having background color
            std::vector<std::vector<Pixel>> pixels(2);
            fmt::println("Image 1: {} pixels", getPixelsNotHavingBackground(image1, colorBackground, pixels[0]));
            fmt::println("Image 2: {} pixels", getPixelsNotHavingBackground(image2, colorBackground, pixels[1]));
            // now we need to make sure the list have the same number of elements
            int index_to_resize = 0;
            size_t resizeToElements = pixels[1].size();
            if (pixels[0].size() > pixels[1].size()) {
                index_to_resize = 1;
                resizeToElements = pixels[0].size();
            }

            std::random_device rd;  // obtain a random seed from hardware
            std::mt19937 gen(rd());  // seed the generator
            std::uniform_int_distribution<> distr(0, pixels[index_to_resize].size() - 1);  // define the range


            fmt::println("Adding {} random pixels for interpolation",
                         std::abs((int) (pixels[0].size() - pixels[1].size())));
            while (pixels[index_to_resize].size() < resizeToElements) {
                size_t copyIndex = distr(gen);
                pixels[index_to_resize].emplace_back(pixels[index_to_resize][copyIndex]);
            }


            // random shuffle of starting positions
            if (shuffleStartingPositions) {
                std::shuffle(std::begin(pixels[0]), std::end(pixels[0]), gen);
            }

            float stepSize = 1.0f / frames;
            float rotateStepSize = 360.0f / frames;

            for (int s = 0; s < frames; s++) {
                std::string fileOut = fileTarget;
                float fraction = 0.0f + (stepSize * s);
                float degrees = 0.0f + (rotateStepSize * s);
                size_t frameNumber = startFrame + s;
                std::string frameNumberString = fmt::format("{:0{}d}", frameNumber, 9);
                fileOut += frameNumberString + ".png";
                fmt::print("generating {}\n", fileOut);
                // Interpolate all pixel colors and pixel positions into a new pixels vector / image
                Image imageTarget(image1.width, image1.height);
                imageTarget.setAlpha(255);
                // check for edge case where we dont need to interpolate
                bool sourceIsTarget = false;
                if (s == 0) {
                    // first frame copy content
                    imageTarget = image1;
                    sourceIsTarget = true;
                } else if (s == frames - 1) {
                    // last frame copy content
                    imageTarget = image2;
                    sourceIsTarget = true;
                } else {
                    // interpolation required
                    for (size_t i = 0; i < pixels[0].size(); i++) {
                        Pixel pixel = interpolate(pixels[0][i], pixels[1][i], fraction);
                        setPixelAdd(imageTarget, pixel.x, pixel.y, pixel.color);
                    }
                    std::vector<Pixel> interpolatedPixels;
                    // store Pixels
                    getPixelsNotHavingBackground(imageTarget, colorBackground, interpolatedPixels);
                    // add motion blur
                    gaussianBlurRGBA(imageTarget);
                    // redraw interpolated Pixels
                    for (size_t i = 0; i < interpolatedPixels.size(); i++) {
                        setPixel(imageTarget, interpolatedPixels[i].x, interpolatedPixels[i].y,
                                 interpolatedPixels[i].color);
                    }
                }
                if (rotate && !sourceIsTarget) {
                    rotateImage(imageTarget, colorBackground, degrees);
                }
                // draw VU meters
                Color colorVUmeter = { 255, 255, 0, 255};
                int vux1 = 20;
                int vux2 = vux1 + (20 * dsp.vumeters[2][frameNumber] );
                int vuy = 20;
                drawLine(imageTarget, vux1, vuy, vux2, vuy, colorVUmeter);
                drawLine(imageTarget, vux1, vuy+1, vux2, vuy+1, colorVUmeter);

                vux2 = vux1 + (20 * dsp.vumeters[1][frameNumber] );
                vuy = 30;
                drawLine(imageTarget, vux1, vuy, vux2, vuy, colorVUmeter);
                drawLine(imageTarget, vux1, vuy+1, vux2, vuy+1, colorVUmeter);

                vux2 = vux1 + (20 * dsp.vumeters[0][frameNumber] );
                vuy = 40;
                drawLine(imageTarget, vux1, vuy, vux2, vuy, colorVUmeter);
                drawLine(imageTarget, vux1, vuy+1, vux2, vuy+1, colorVUmeter);
                fmt::println("VU: {} {} {}", dsp.vumeters[2][frameNumber], dsp.vumeters[1][frameNumber], dsp.vumeters[0][frameNumber]);

                imageTarget.savePNG(fileOut);
            }

        } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }

        return 0;
    }
}

int main() {

    float fc1 = 0.1f; // Low-pass filter cutoff frequency
    float fl = 0.1f; // Mid-pass filter lower cutoff frequency
    float fh = 0.5f; // Mid-pass filter higher cutoff frequency
    float fc2 = 0.5f; // High-pass filter cutoff frequency
    mvm::Dsp dsp;
    if (dsp.loadOgg("shona.ogg")) {
        dsp.generateDefaultBands(fc1, fl, fh, fc2);
        dsp.generateVuMeters(25);
    }

    std::string img1 = mvm::walker(0.56);
    std::string img2 = mvm::walker(0.596);
    std::string imgInterpolated = "walker2_anim_";
    int frames = 100;
    int startFrame = 1;
    mvm::Color colorBackground = {0, 0, 0, 255};
    mvm::interpolate(dsp,img1, img2, imgInterpolated, frames, startFrame, colorBackground, false);

    startFrame += frames;
    img1 = img2;
    img2 = mvm::walker(0.58);
    mvm::interpolate(dsp,img1, img2, imgInterpolated, frames, startFrame, colorBackground, false);

    startFrame += frames;
    img1 = img2;
    img2 = mvm::walker(0.59);
    mvm::interpolate(dsp,img1, img2, imgInterpolated, frames, startFrame, colorBackground, true);

    startFrame += frames;
    img1 = img2;
    img2 = mvm::walker(0.66);
    mvm::interpolate(dsp,img1, img2, imgInterpolated, frames, startFrame, colorBackground, false);

    startFrame += frames;
    img1 = img2;
    img2 = mvm::walker(0.33);
    mvm::interpolate(dsp,img1, img2, imgInterpolated, frames, startFrame, colorBackground, true);

    startFrame += frames;
    img1 = img2;
    img2 = mvm::walker(0.22);
    mvm::interpolate(dsp,img1, img2, imgInterpolated, frames, startFrame, colorBackground, false);

    startFrame += frames;
    img1 = img2;
    img2 = mvm::walker(0.57);
    mvm::interpolate(dsp,img1, img2, imgInterpolated, frames, startFrame, colorBackground, false);

    startFrame += frames;
    img1 = img2;
    img2 = mvm::walker(0.56);
    mvm::interpolate(dsp,img1, img2, imgInterpolated, frames, startFrame, colorBackground, true);

    return 0;
}
