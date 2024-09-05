#include <drogon/drogon.h>
#include "../databases/comparator.hpp"
int main() {
    //Set HTTP listener address and port
    //drogon::app().addListener("0.0.0.0",8081);
    //Load config file
    drogon::app().loadConfigFile(getProjectRootPath() + "/config.json");
    drogon::app().setClientMaxBodySize(10737418240); // 10GB
    //Run HTTP framework,the method will block in the internal event loop
    drogon::app().run();
    return 0;
}
