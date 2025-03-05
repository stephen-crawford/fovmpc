#ifndef CATMULLROMSPLINE_H
#define CATMULLROMSPLINE_H

#include "CatmullRomSegment.h"
#include <vector>
#include <stdexcept>


// https://www.desmos.com/calculator/9kazaxavsf

// https://www.desmos.com/calculator/cahqdxeshd
namespace splines {

    template <typename T, unsigned int DIM>
    class CatmullRomSpline {
    public:
        using VectorDIM = typename CatmullRomSegment<T, DIM>::VectorDIM;

        // Constructor from a list of segments
        explicit CatmullRomSpline(const std::vector<CatmullRomSegment<T, DIM>>& segments);

        // Evaluate the spline at a given parameter t
        VectorDIM eval(T t) const;

        // Append a new segment
        void appendSegment(const CatmullRomSegment<T, DIM>& segment);

        // Get number of segments
        std::size_t numSegments() const;

        // Get all segments
        const std::vector<CatmullRomSegment<T, DIM>>& getSegments() const;

    private:
        std::vector<CatmullRomSegment<T, DIM>> segments_;
    };

    // Implementation must be included in the header for templates
    template<typename T, unsigned int DIM>
    CatmullRomSpline<T, DIM>::CatmullRomSpline(const std::vector<CatmullRomSegment<T, DIM>>& segments)
        : segments_(segments) {
        if (segments_.empty()) {
            throw std::runtime_error("CatmullRomSpline requires at least one segment.");
        }
    }

    template<typename T, unsigned int DIM>
    typename CatmullRomSpline<T, DIM>::VectorDIM
    CatmullRomSpline<T, DIM>::eval(T t) const {
        if (t < 0 || t > static_cast<T>(segments_.size())) {
            throw std::runtime_error("CatmullRomSpline: t is out of range.");
        }

        size_t segment_index = static_cast<size_t>(t);
        if (segment_index >= segments_.size()) {
            segment_index = segments_.size() - 1;
        }
        T local_t = t - static_cast<T>(segment_index);

        return segments_[segment_index].eval(local_t);
    }

    template<typename T, unsigned int DIM>
    void CatmullRomSpline<T, DIM>::appendSegment(const CatmullRomSegment<T, DIM>& segment) {
        segments_.push_back(segment);
    }

    template<typename T, unsigned int DIM>
    size_t CatmullRomSpline<T, DIM>::numSegments() const {
        return segments_.size();
    }

    template<typename T, unsigned int DIM>
    const std::vector<CatmullRomSegment<T, DIM>>&
    CatmullRomSpline<T, DIM>::getSegments() const {
        return segments_;
    }

} // namespace splines

#endif // CATMULLROMSPLINE_H
