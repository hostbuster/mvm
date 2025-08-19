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

#include <vorbis/vorbisfile.h>

bool loadOgg(const std::string& filename, std::vector<float>& pcmData, int& channels, int& sampleRate)
{
    // Open the file
    FILE* file = fopen(filename.c_str(), "rb");
    if (!file)
    {
        std::cerr << "Failed to open file " << filename << std::endl;
        return false;
    }

    // Initialize the Ogg Vorbis stream
    OggVorbis_File oggFile;
    if (ov_open(file, &oggFile, nullptr, 0) < 0)
    {
        std::cerr << "Failed to open Ogg file " << filename << std::endl;
        fclose(file);
        return false;
    }

    // Get information about the audio stream
    vorbis_info* info = ov_info(&oggFile, -1);
    channels = info->channels;
    sampleRate = info->rate;

    // Decode the audio data and store it in a vector
    long totalSamples = ov_pcm_total(&oggFile, -1);
    pcmData.resize(totalSamples * channels);
    float** pcm = new float*[channels];
    for (int i = 0; i < channels; i++)
    {
        pcm[i] = &pcmData[i * totalSamples];
    }
    long samplesRead = 0;
    int bitstream = 0;
    while (samplesRead < totalSamples)
    {
        long samplesToRead = std::min(totalSamples - samplesRead, 4096L);
        long bytesRead = ov_read_float(&oggFile, pcm, samplesToRead, &bitstream);
        if (bytesRead <= 0)
        {
            break;
        }
        samplesRead += bytesRead;
        pcm += bytesRead * channels;
    }
    delete[] pcm;

    // Clean up and return success
    ov_clear(&oggFile);
    fclose(file);
    return true;
}



int main()
{
    std::string filename = "example.ogg";
    std::vector<float> pcmData;
    int channels, sampleRate;
    if (loadOgg(filename, pcmData, channels, sampleRate))
    {
        std::cout << "Loaded Ogg file " << filename << " with " << channels << " channels and sample rate " << sampleRate << std::endl;
    }
    else
    {
        std::cerr << "Failed to load Ogg file " << filename << std::endl;
    }


    // filters
    std::vector<float> input = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f };
    std::vector<float> output;
    float fc1 = 0.1f; // Low-pass filter cutoff frequency
    float fl = 0.1f; // Mid-pass filter lower cutoff frequency
    float fh = 0.5f; // Mid-pass filter higher cutoff frequency
    float fc2 = 0.5f; // High-pass filter cutoff frequency

    lowpassFilter(input, output, fc1);
// output: { 1.0f, 1.09093f, 1.25968f, 1.49901f, 1.79858f }

    midpassFilter(input, output, fl, fh);
// output: { 0.0f, 0.904508f, 1.788854f, 2.671572f, 3.536355f }

    highpassFilter(input, output, fc2);
// output: { 0.0f, 0.0956614f,

    return 0;
}


#include <vector>

// Low-pass filter with cutoff frequency fc
void lowpassFilter(const std::vector<float>& input, std::vector<float>& output, float fc)
{
    float c = 1.0f / (1.0f + 2.0f * 3.14159265f * fc);
    output.resize(input.size());
    output[0] = input[0];
    for (int i = 1; i < input.size(); i++)
    {
        output[i] = c * (input[i] + input[i - 1]) + (1.0f - 2.0f * c) * output[i - 1];
    }
}

// Mid-pass filter with cutoff frequencies fl and fh
void midpassFilter(const std::vector<float>& input, std::vector<float>& output, float fl, float fh)
{
    float c1 = 2.0f * 3.14159265f * fl;
    float c2 = 2.0f * 3.14159265f * fh;
    float c = 1.0f / (1.0f + c1 / c2 + c1 * c1 / (c2 * c2));
    output.resize(input.size());
    output[0] = input[0];
    for (int i = 1; i < input.size(); i++)
    {
        output[i] = c * (input[i] - input[i - 1] + c1 / c2 * output[i - 1]) + (1.0f - c * (1.0f + c1 / c2)) * output[i - 1];
    }
}

// High-pass filter with cutoff frequency fc
void highpassFilter(const std::vector<float>& input, std::vector<float>& output, float fc)
{
    float c = 1.0f / (1.0f + 2.0f * 3.14159265f * fc);
    output.resize(input.size());
    output[0] = input[0];
    for (int i = 1; i < input.size(); i++)
    {
        output[i] = c * (input[i] - input[i - 1]) + (1.0f - c) * output[i - 1];
    }
}

#include <vector>

// Band-pass filter with center frequency fc and bandwidth bw
void bandpassFilter(const std::vector<float>& input, std::vector<float>& output, float fc, float bw)
{
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
    for (int i = 2; i < input.size(); i++)
    {
        output[i] = c1 * input[i] + c2 * input[i - 1] + c3 * input[i - 2];
    }
}

std::vector<double> lowpassFilter(const std::vector<double>& input, double cutoff, double sampleRate, double resonance)
{
    double c = 1.0 / tan(M_PI * cutoff / sampleRate);
    double a1 = 1.0 / (1.0 + resonance * c + c * c);
    double a2 = 2 * a1;
    double a3 = a1;
    double b1 = 2.0 * (1.0 - c*c) * a1;
    double b2 = (1.0 - resonance * c + c*c) * a1;

    std::vector<double> output(input.size());

    // Direct Form I implementation
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    for (int i = 0; i < input.size(); i++) {
        output[i] = a1 * input[i] + a2 * x1 + a3 * x2 - b1 * y1 - b2 * y2;
        x2 = x1;
        x1 = input[i];
        y2 = y1;
        y1 = output[i];
    }

    return output;
}








// Apply a series of band-pass filters with different center frequencies and bandwidths
void equalize(const std::vector<float>& input, std::vector<float>& output)
{
    output.resize(input.size());
    std::vector<float> temp(input.size());
    float fc[] = { 50.0f, 100.0f, 200.0f, 400.0f, 800.0f };
    float bw[] = { 20.0f, 40.0f, 80.0f, 160.0f, 320.0f };
    for (int i = 0; i < 5; i++)
    {
        bandpassFilter(input, temp, fc[i], bw[i]);
        for (int j = 0; j < input.size(); j++)
        {
            output[j] += temp[j];
        }
    }
}

