#include "api_v1_Server.h"

using namespace api::v1;
using namespace std;

std::vector<std::pair<uint32_t, int32_t>> parseStringUI(const std::string &s)
{
    std::vector<std::pair<uint32_t, int32_t>> result;
    std::stringstream ss(s);
    char c1, c2, comma;
    int x, y;

    ss >> c1; // read [
    while (ss >> c1 >> x >> comma >> y >> c2)
    {
        if (c1 != '(' || comma != ',' || c2 != ')')
        {
            std::cerr << "Invalid format: " << s << std::endl;
            return {};
        }
        result.emplace_back(x, y);

        ss >> c1;
        if (c1 == ']')
        {
            break;
        }
    }
    if (c1 != ']')
    {
        std::cerr << "Invalid format: " << s << std::endl;
        return {};
    }
    return result;
}

std::vector<std::pair<uint32_t, uint32_t>> parseStringUU(const std::string &s)
{
    std::vector<std::pair<uint32_t, uint32_t>> result;
    std::stringstream ss(s);
    char c1, c2, comma;
    int x, y;

    ss >> c1; // read [
    while (ss >> c1 >> x >> comma >> y >> c2)
    {
        if (c1 != '(' || comma != ',' || c2 != ')')
        {
            std::cerr << "Invalid format: " << s << std::endl;
            return {};
        }
        result.emplace_back(x, y);

        ss >> c1;
        if (c1 == ']')
        {
            break;
        }
    }
    if (c1 != ']')
    {
        std::cerr << "Invalid format: " << s << std::endl;
        return {};
    }
    return result;
}

void print_vector(const std::vector<std::string> &v)
{
    for (const std::string &s : v)
    {
        std::cout << s << std::endl;
    }
}

Server::Server() : squid(constants::SmallFastComp, true)
{
    api_keys = std::unordered_set<std::string>{MasterApiKey};
    std::string path = getProjectRootPath();
    squid.setVCFData(path + "/data/sample.vcf");
    squid.setPhenoData(path + "/data/sample.csv");

    squid.printMeta();

    squid.printData(true);
    squid.printPhenoData(true);

    std::cout << "Server is ready" << std::endl;
};

void Server::getContext(const HttpRequestPtr &req,
                        std::function<void(const HttpResponsePtr &)> &&callback,
                        const std::string &apikey) const
{
    LOG_DEBUG << "Running Get Context query with from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::ostringstream stream;
    helib::Context &copy = squid.getMeta().data->context;

    copy.writeToJSON(stream);
    ret["result"] = stream.str();
    auto resp = HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
    return;
}

void Server::printDB(const HttpRequestPtr &req,
                     std::function<void(const HttpResponsePtr &)> &&callback,
                     const std::string &apikey) const
{
    LOG_DEBUG << "Running Print DB query with from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    ret["result"] = squid.printDBString(true);
    auto resp = HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
    return;
}

void Server::authorizeAPI(const HttpRequestPtr &req,
                          std::function<void(const HttpResponsePtr &)> &&callback,
                          const std::string &apikey)
{
    LOG_DEBUG << "Running authorize API from a user with Key:" << apikey;

    std::cout << "Running authorize API from a user with Key:" << apikey << std::endl;

    Json::Value ret;

    if (api_keys.count(apikey) == 0)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::string temp_file = "temp";
    HttpRequest *req_ptr = req.get();
    std::ofstream outfile(temp_file, std::ios::out);
    outfile.write(req_ptr->bodyData(), req_ptr->bodyLength());
    outfile.close();

    std::ifstream in_pubkey_file;
    in_pubkey_file.open(temp_file, std::ios::in);
    if (in_pubkey_file.is_open())
    {
        helib::PubKey client_public_key = helib::PubKey::readFromJSON(in_pubkey_file, squid.getContext());
        in_pubkey_file.close();
        addKSK(client_public_key, apikey);
        std::cout << "Added key to store" << std::endl;
        ret["result"] = "success";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }
    else
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }
}

void Server::countingQueryAPI(const HttpRequestPtr &req,
                              std::function<void(const HttpResponsePtr &)> &&callback,
                              std::string query,
                              std::string conj,
                              const std::string &apikey) const
{
    LOG_DEBUG << "Running counting query with " << query << " from user with API Key: " << apikey;

    Json::Value ret;

    if (apikey != MasterApiKey)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::vector<std::pair<uint32_t, uint32_t>> pairs = parseStringUU(query);

    if (pairs.size() == 0)
    {
        ret["result"] = "query failed to parse or is empty";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    bool conjunctive;

    if (conj == "0")
    {
        conjunctive = false;
    }
    else if (conj == "1")
    {
        conjunctive = true;
    }
    else
    {
        ret["result"] = "conjunctive failed to parse or is empty";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    helib::Ctxt result = squid.countQuery(conjunctive, pairs);

    auto ksk = key_switch_store.at(apikey);

    result.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));

    std::stringstream ss;
    result.writeToJSON(ss);
    ret["result"] = ss.str();

    auto resp = HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
}

void Server::mafQueryAPI(const HttpRequestPtr &req,
                         std::function<void(const HttpResponsePtr &)> &&callback,
                         std::string query,
                         std::string conj,
                         std::string target,
                         const std::string &apikey) const
{
    LOG_DEBUG << "Running MAF query with " << query << " from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::vector<std::pair<uint32_t, uint32_t>> pairs = parseStringUU(query);

    if (pairs.size() == 0)
    {
        ret["result"] = "query failed to parse or is empty";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    bool conjunctive;

    if (conj == "0")
    {
        conjunctive = false;
    }
    else if (conj == "1")
    {
        conjunctive = true;
    }
    else
    {
        ret["result"] = "conjunctive failed to parse or is empty";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    int i_target;
    try
    {
        i_target = std::stoi(target);
    }
    catch (const std::invalid_argument &e)
    {
        std::cout << "Invalid argument: " << e.what() << std::endl;

        ret["result"] = "couldn't convert target to int";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }
    catch (const std::out_of_range &e)
    {
        std::cout << "Out of range: " << e.what() << std::endl;
        ret["result"] = "target int was out of range";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    helib::Ctxt result = squid.MAFQuery(i_target, conjunctive, pairs);

    auto ksk = key_switch_store.at(apikey);

    result.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));

    std::stringstream ss;
    result.writeToJSON(ss);
    ret["result"] = ss.str();

    auto resp = HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
}

void Server::PRSQueryAPI(const HttpRequestPtr &req,
                         std::function<void(const HttpResponsePtr &)> &&callback,
                         std::string params,
                         const std::string &apikey) const
{
    LOG_DEBUG << "Running PRS query with " << params << " from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::vector<std::pair<uint32_t, int32_t>> pairs = parseStringUI(params);

    if (pairs.size() == 0)
    {
        ret["result"] = "query failed to parse or is empty";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::vector<helib::Ctxt> results = squid.PRSQuery(pairs);
    uint32_t num_rows = squid.getNumRows();

    auto ksk = key_switch_store.at(apikey);

    for (size_t i = 0; i < results.size(); i++)
    {
        results[i].PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));
    }
    for (size_t i = 0; i < results.size(); i++)
    {
        std::stringstream ss;
        results[i].writeToJSON(ss);
        ret["result_" + std::to_string(i)] = ss.str();
    }
    ret["num_rows"] = num_rows;
    ret["num_cipher_texts"] = results.size();

    auto resp = HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
}

void Server::similarityQueryAPI(const HttpRequestPtr &req,
                                std::function<void(const HttpResponsePtr &)> &&callback,
                                std::string threshold,
                                std::string binaryPheno,
                                const std::string &apikey) const

{
    LOG_DEBUG << "Running Similarity query with " << threshold << " threshold and " << binaryPheno << " phenotype column from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::string temp_file = "temp";
    HttpRequest *req_ptr = req.get();
    std::ofstream outfile(temp_file, std::ios::out);
    outfile.write(req_ptr->bodyData(), req_ptr->bodyLength());
    outfile.close();

    std::ifstream in_ctxt_stream;

    vector<helib::Ctxt> ctxts;
    in_ctxt_stream.open(temp_file, std::ios::in);
    if (in_ctxt_stream.is_open())
    {
        for (int i = 0; i < squid.num_snp_cols; i++)
        {
            helib::Ctxt ctxt = helib::Ctxt::readFromJSON(in_ctxt_stream, squid.getMeta().data->publicKey);
            ctxts.push_back(ctxt);
        }
        in_ctxt_stream.close();
    }
    else
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    pair<helib::Ctxt, helib::Ctxt> result = squid.similarityQuery(std::stoi(binaryPheno), ctxts, std::stoi(threshold));
    
    auto ksk = key_switch_store.at(apikey);
    result.first.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));
    result.second.PublicKeySwitch(std::make_pair(std::ref(ksk.first), std::ref(ksk.second)));

    std::stringstream ss;
    result.first.writeToJSON(ss);
    ret["result_with"] = ss.str();

    std::stringstream ss2;
    result.second.writeToJSON(ss2);
    ret["result_without"] = ss2.str();

    ret["result"] = "success";
    auto resp = HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
    return;
}

void Server::getHeadersAPI(const HttpRequestPtr &req,
                           std::function<void(const HttpResponsePtr &)> &&callback,
                           const std::string &apikey) const
{
    LOG_DEBUG << "Running get headers with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    string s = "SNPs: {";

    vector<string> column_headers = squid.getColumnHeaders();

    for (int i = 0; i < column_headers.size() - 1; i++)
    {
        s += to_string(i) + ":" + column_headers[i] + ", ";
    }
    s += to_string(column_headers.size() - 1) + ":" + column_headers[column_headers.size() - 1] + "}";

    s += "\n Binary Phenotypes: {";

    vector<string> binary_pheno_headers = squid.getBinaryPhenoHeaders();

    for (int i = 0; i < binary_pheno_headers.size() - 1; i++)
    {
        s += to_string(i) + ":" + binary_pheno_headers[i] + ", ";
    }
    s += to_string(binary_pheno_headers.size() - 1) + ":" + binary_pheno_headers[binary_pheno_headers.size() - 1] + "}";

    s += "\n Continuous Phenotypes: {";

    vector<string> continuous_pheno_headers = squid.getContinuousPhenoHeaders();

    for (int i = 0; i < continuous_pheno_headers.size() - 1; i++)
    {
        s += to_string(i) + ":" + continuous_pheno_headers[i] + ", ";
    }
    s += to_string(continuous_pheno_headers.size() - 1) + ":" + continuous_pheno_headers[continuous_pheno_headers.size() - 1] + "}";

    ret["result"] = s;
    auto resp = HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
}

void Server::getOwnerPublicKey(const HttpRequestPtr &req,
                               std::function<void(const HttpResponsePtr &)> &&callback,
                               const std::string &apikey) const
{
    LOG_DEBUG << "Running Get Owner Public Key with from user with API Key: " << apikey;

    Json::Value ret;

    if (api_keys.count(apikey) == 0)
    {
        ret["result"] = "failed";
        auto resp = HttpResponse::newHttpJsonResponse(ret);
        callback(resp);
        return;
    }

    std::ostringstream out;

    squid.getMeta().data->publicKey.writeToJSON(out);

    ret["result"] = out.str();

    auto resp = HttpResponse::newHttpJsonResponse(ret);
    callback(resp);
    return;
}

void Server::addKSK(helib::PubKey &client_pk, string id)
{
    auto ksk = client_pk.genPublicKeySwitchingKey(squid.getSecretKey());
    key_switch_store.emplace(id, ksk);
}