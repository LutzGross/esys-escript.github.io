
/*****************************************************************************
*
* Copyright (c) 2003-2017 by The University of Queensland
* http://www.uq.edu.au
*
* Primary Business: Queensland, Australia
* Licensed under the Apache License, version 2.0
* http://www.apache.org/licenses/LICENSE-2.0
*
* Development until 2012 by Earth Systems Science Computational Center (ESSCC)
* Development 2012-2013 by School of Earth Sciences
* Development from 2014 by Centre for Geoscience Computing (GeoComp)
*
*****************************************************************************/

#ifndef __ESCRIPT_PARAMS_H__
#define __ESCRIPT_PARAMS_H__

#include <boost/python/list.hpp>

#include <string>
#include <unordered_set>

namespace escript
{

class EscriptParams
{
    typedef std::unordered_set<std::string> FeatureSet;

public:
    EscriptParams();

    int getInt(const std::string& name, int sentinel = 0) const;
    void setInt(const std::string& name, int value);
    boost::python::list listEscriptParams() const;

    inline int getAutoLazy() const { return autoLazy; }
    inline int getLazyStrFmt() const { return lazyStrFmt; }
    inline int getLazyVerbose() const { return lazyVerbose; }
    inline int getResolveCollective() const { return resolveCollective; }
    inline int getTooManyLevels() const { return tooManyLevels; }
    inline int getTooManyLines() const { return tooManyLines; }

    bool hasFeature(const std::string& name) const;
    boost::python::list listFeatures() const;

private:
    FeatureSet features;
    // the number of parameters is small enough to avoid a map for performance
    // reasons
    int autoLazy;
    int lazyStrFmt;
    int lazyVerbose;
    int resolveCollective;
    int tooManyLevels;
    int tooManyLines;
};


extern EscriptParams escriptParams;

/**
    \brief Set the value of a named parameter.
    See listEscriptParams() for available parameters.
*/
inline void setEscriptParamInt(const std::string& name, int value)
{
   escriptParams.setInt(name, value);
}

/**
    \brief get the value of a named parameter.
    See listEscriptParams() for available parameters.
*/
inline int getEscriptParamInt(const std::string& name, int sentinel=0)
{
    return escriptParams.getInt(name, sentinel);
}

/**
    \brief describe available parameters.
    \return a list of tuples (parameter name, value, description)
*/
inline boost::python::list listEscriptParams()
{
    return escriptParams.listEscriptParams();
}

/**
    \brief returns true if escript was compiled with the feature `name`,
           false otherwise.
*/
inline bool hasFeature(const std::string& name)
{
    return escriptParams.hasFeature(name);
}

/**
    \brief returns a list of features escript was compiled with.
    \return a boost python list of strings
*/
inline boost::python::list listFeatures()
{
    return escriptParams.listFeatures();
}

} // namespace escript

#endif // __ESCRIPT_PARAMS_H__

