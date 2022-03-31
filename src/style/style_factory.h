#pragma once

#include <map>
#include <string>
#include <iostream>

template<typename T, typename ...Args>
class StyleFactory:
  private std::map<std::string, T *(*)(Args ...)>
{
public:
  template<typename U>
  void
  register_class(const std::string &key)
  {
    try_emplace(key, [](Args ...args) -> T *{ return new U(args...); });
  }

  T *
  new_instance(const std::string &key, Args ...args)
  {
    const const_iterator &constructor = find(key);

    if (constructor == cend())
    {
      std::cout << "Error: key \"" << key << "\" does not exist."  << std::endl;
      return nullptr;
    }

    return constructor->second(args...);
  }
};
