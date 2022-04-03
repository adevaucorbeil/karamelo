#pragma once

#include <map>
#include <string>
#include <iostream>

template<typename T, typename ...Args>
class StyleFactory
{
  std::map<std::string, T *(*)(Args ...)> internal_map;

public:
  template<typename U>
  void
  register_class(const std::string &key)
  {
    internal_map[key] = [](Args ...args) -> T *{ return new U(args...); };
  }

  T *
  new_instance(const std::string &key, Args ...args)
  {
    const typename std::map<std::string, T *(*)(Args ...)>::const_iterator &constructor = internal_map.find(key);

    if (constructor == internal_map.cend())
    {
      std::cout << "Error: key \"" << key << "\" does not exist." << std::endl;
      return nullptr;
    }

    return constructor->second(args...);
  }
};
