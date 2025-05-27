#include <bits/stdc++.h>

using namespace std;

class Node;

class Element;

class Node {
private:
    string name;
    double voltage;
    vector<Element *> connected_elements;
public:
    Node(string name_) : name(name_) {}

    vector<Element *> getConnectedElements() { return connected_elements; }

    string getName() { return name; }

    void addConnectedElement(Element *e) {
        connected_elements.push_back(e);
    }

    void removeConnectedElement(Element *e) {
        auto x = find(connected_elements.begin(), connected_elements.end(), e);
        if (x != connected_elements.end())
            connected_elements.erase(x);
    }

    void setName(string new_name) {
        name = new_name;
    }
};

class Circuit {
protected:
    vector<Element *> elements;
    vector<Node *> nodes;
    map<string, Node *> node_access;
    Node *gnd=nullptr;
public:
    vector<Element *> &getElements() { return elements; }

    map<string, Node *> getNodeAccess() { return node_access; }

    Node *getCreateNode(string name) {
        if (!node_access.count(name)) {
            Node *newnode = new Node(name);
            nodes.push_back(newnode);
            node_access[name] = newnode;
        }
        return node_access[name];
    }

    void addElement(Element *element) {
        elements.push_back(element);
    }

    Node* getGround() {return gnd;}

    bool setGround(string name){
        if(!node_access.count(name))
            return false;
        gnd=node_access[name];
        return true;
    }

    bool removeGround(string name){
        if(gnd!=node_access[name])
            return false;
        gnd=nullptr;
        return true;
    }

    bool renameNode(string old_name, string new_name, string& err) {
        if (node_access.find(old_name) == node_access.end()) {
            err = "ERROR: Node " + old_name + " does not exist in the circuit\n";
            return false;
        }
        if (node_access.find(new_name) != node_access.end()) {
            err = "ERROR: Node name " + new_name + " already exists\n";
            return false;
        }
        Node* node = node_access[old_name];
        node->setName(new_name);
        node_access.erase(old_name);
        node_access[new_name] = node;
        err = "SUCCESS: Node renamed from " + old_name + " to " + new_name + "\n";
        return true;
    }

};

class Element {
protected:
    string name;
    double value;
    Node *node1;
    Node *node2;
public:
    Element(string name_, double value_, Node *n1, Node *n2) :
            name(name_), value(value_), node1(n1), node2(n2) {}

    string getName() {
        return name;
    }

    virtual string getType() = 0;

    double getValue() { return value; }

    Node *getFirstNode() { return node1; }

    Node *getSecondNode() { return node2; }

    virtual ~Element() = default;
};

class Resistor : public Element {
public:
    Resistor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Resistor"; }
};

class Capacitor : public Element {
public:
    Capacitor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Capacitor"; }
};

class Inductor : public Element {
public:
    Inductor(string name_, double value_, Node *n1, Node *n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override { return "Inductor"; }
};

class Diode : public Element {
protected:
    string model;
public:
    Diode(string name, Node *n1, Node *n2) :
            Element(name, 0, n1, n2), model("D") {}

    string getType() override { return "Diode"; }

    virtual double calculateCurrent(double voltage) {

    }
};

class ZenerDiode : public Diode {
private:
    double breakdownvoltage;
public:
    ZenerDiode(string name, Node *n1, Node *n2, double breakdownvoltage_ = 5)
            : Diode(name, n1, n2), breakdownvoltage(breakdownvoltage_) {
        model = "Z";
    }

    string getType() override { return "Zener Duiode"; }

    double calculateCurrent(double voltage) override {

    }
};

class Controller {
public:
    string handleError(string name, Circuit *circuit, double value) {
        string type;
        switch (name[0]) {
            case 'R':
                type = "Resistor";
                break;
            case 'C':
                type = "Capacitor";
                break;
            case 'L':
                type = "Inductor";
                break;
            case 'D':
                type = "Diode";
                break;
        }
        if (value <= 0 && (type == "Resistor" || type == "Capacitor" || type == "Inductor"))
            return "Error: " + type + " cannot be zero or negative\n";
        for (auto element: circuit->getElements()) {
            if (!element)
                continue;
            if (element->getName() == name) {
                switch (name[0]) {
                    case 'R':
                        return "Resistor " + name + " already exists in the circuit\n";
                    case 'C':
                        return "Capacitor " + name + " already exists in the circuit\n";
                    case 'L':
                        return "Inductor " + name + " already exists in the circuit\n";
                    case 'D':
                        return "Diode" + name + " already exists in the circuit\n";
                }
            }
        }
        return "";
    }

    string addNewElement(string node1, string node2, string name, double value, Circuit *circuit) {
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        Element *element = nullptr;
        switch (name[0]) {
            case 'R':
                element = new Resistor(name, value, n1, n2);
                break;
            case 'L':
                element = new Inductor(name, value, n1, n2);
                break;
            case 'C':
                element = new Capacitor(name, value, n1, n2);
                break;
        }
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }

    string removeElement(string name, Circuit *circuit) {
        auto &elements = circuit->getElements();
        for (auto it = elements.begin(); it != elements.end(); it++) {
            if ((*it)->getName() == name) {
                (*it)->getFirstNode()->removeConnectedElement(*it);
                (*it)->getSecondNode()->removeConnectedElement(*it);
                Element *deleting = *it;
                elements.erase(it);
                delete deleting;
                return name + " removed successfully!\n";
            }
        }
        switch (name[0]) {
            case 'R':
                return "Error: Cannot delete resistor; component not found\n";
            case 'L':
                return "Error: Cannot delete inductor; component not found\n";
            case 'C':
                return "Error: Cannot delete capacitor; component not found\n";
            case 'D':
                return "Error: Cannot delete diode; component not found\n";
        }
    }

    string addDiode(string node1, string node2, string name, string model, Circuit *circuit) {
        Node *n1 = circuit->getCreateNode(node1);
        Node *n2 = circuit->getCreateNode(node2);
        if (model == "Z") {
            ZenerDiode *zenerDiode = new ZenerDiode(name, n1, n2);
            circuit->addElement(zenerDiode);
            n1->addConnectedElement(zenerDiode);
            n2->addConnectedElement(zenerDiode);
            return "Zener Diode added successfully!\n";
        } else if (model == "D") {
            Diode *diode = new Diode(name, n1, n2);
            circuit->addElement(diode);
            n1->addConnectedElement(diode);
            n2->addConnectedElement(diode);
            return "Diode added successfully!\n";
        }
    }

    string addGround(string node, Circuit* circuit){
        if(circuit->getGround())
            return "Error! ground node already exists in the circuit\n";
        Node *n1 = circuit->getCreateNode(node);
        circuit->setGround(node);
        return "Ground added to node " + node +" successfully!\n";
    }

    string deleteGround(string node, Circuit* circuit){
        if(!circuit->getNodeAccess().count(node))
            return "Node does not exist\n";
        if(!circuit->removeGround(node))
            return "Error! Node is not set as ground\n";
        return "Ground removed from node "+node+" successfully!\n";
    }

    void showNodes(Circuit* circuit){
        cout<<"Available nodes:\n";
        if(!circuit->getNodeAccess().size())
            cout<<"No node exists in the circuit\n";
        int index=1;
        for(auto node: circuit->getNodeAccess()){
            cout<<node.first;
            if(index!=circuit->getNodeAccess().size())
                cout<<", ";
            else
                cout<<"\n";
            index++;
        }
    }

    void list(Circuit* circuit){
        cout<<"Elements:\n";
        if(!circuit->getElements().size())
            cout<<"No elements exists in the circuit\n";
        for(auto element: circuit->getElements()){
            cout<<element->getType()<<": "<<element->getName()<<", value: "<<element->getValue()<<endl;
        }
    }

    void listElement(string type, Circuit* circuit){
        switch(type[0]){
            case 'R':
                type = "Resistor";
                break;
            case 'C':
                type = "Capacitor";
                break;
            case 'L':
                type = "Inductor";
                break;
            case 'D':
                type = "Diode";
                break;
        }
        cout<<type<<"(s):\n";
        int index=0;
        for(auto element: circuit->getElements()){
            if(element->getType()==type) {
                cout << element->getType() << ": " << element->getName() << ", value: " << element->getValue() << endl;
                index++;
            }
        }if(!index)
            cout<<"No "<<type<<" exist in the circuit\n";
    }

    void showCircuitDetails(Circuit *circuit) {
        cout << "Circuit Details:\nElements:\n";
        for (auto element: circuit->getElements()) {
            cout << "type: " << element->getType() << "   name: " << element->getName() << "  value: "
                 << element->getValue()
                 << "  node1: " << element->getFirstNode()->getName() << "  node2: "
                 << element->getSecondNode()->getName() << endl;
        }
        cout << "Node Details:\n";
        for (auto node: circuit->getNodeAccess()) {
            cout << node.first << " : " << endl;
            if (node.second->getConnectedElements().size() == 0)
                cout << "No element is connected to this node!\n";
            else {
                for (auto element: node.second->getConnectedElements())
                    cout << element->getName() << endl;
            }
        }
        cout<<"GND:\n";
        if(circuit->getGround())
            cout<<circuit->getGround()->getName()<<endl;
        else
            cout<<"No ground is specified for this circuit\n";
    }

};

class View {
private:
    Circuit *circuit;
    Controller controller;
public:
    void run() {
        circuit = new Circuit();
        string input;
        smatch match;
        regex add_element(R"(^\s*add\s+([A-Za-z])(\w+)\s+(\S+)\s+(\S+)\s+(-?[\d\.]+(?:[eE][+-]?\d+)?)([GMkmunp]?)\s*$)");
        regex remove_element(R"(^\s*delete\s+([A-Za-z])(\w+)\s*$)");
        regex add_diode(R"(^\s*add\s+(D)(\w+)\s+(\S+)\s+(\S+)\s+(D|Z)\s*$)");
        regex show_details(R"(^\s*show\s+details\s*$)");
        regex add_ground(R"(^\s*add\s+(\w+)\s+(\w+)\s*$)");
        regex delete_ground(R"(^\s*delete\s+(\w+)\s+(\w+)\s*$)");
        regex show_nodes(R"(^\s*.nodes\s*$)");
        regex list(R"(^\s*.list\s*$)");
        regex list_element(R"(^\s*.list\s+(\w+)\s*$)");
        regex rename_node(R"(^\s*\.rename\s+node\s+(\S+)\s+(\S+)\s*$)");
        regex exit(R"(^exit$)");
        while (true) {
            getline(cin, input);
            if (regex_match(input, match, add_element)) {
                if (match[1] != "R" && match[1] != "L" && match[1] != "C") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                string name = match[1].str() + match[2].str();
                string node1 = match[3].str();
                string node2 = match[4].str();
                string number = match[5].str();
                string unit = match[6].str();
                double value;
                try {
                    value = stod(number);
                } catch (const invalid_argument &e) {
                    cout << "Error: Invalid numeric value\n";
                    continue;
                }
                string error = controller.handleError(name, circuit, value);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if (!unit.empty()) {
                    switch (unit[0]) {
                        case 'G':
                            value *= 1e9;
                            break;
                        case 'M':
                            value *= 1e6;
                            break;
                        case 'k':
                        case 'K':
                            value *= 1e3;
                            break;
                        case 'm':
                            value *= 1e-3;
                            break;
                        case 'u':
                            value *= 1e-6;
                            break;
                        case 'n':
                            value *= 1e-9;
                            break;
                    }
                }
                cout << controller.addNewElement(node1, node2, name, value, circuit);

            } else if (regex_match(input, match, remove_element)) {
                if (match[1] != "R" && match[1] != "L" && match[1] != "C" && match[1]!="D") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout << controller.removeElement(match[1].str() + match[2].str(), circuit);
            } else if (regex_match(input, match, add_diode)) {
                string name = match[1].str() + match[2].str();
                string node1= match[3];
                string node2=match[4];
                string model=match[5];
                string error = controller.handleError(name, circuit, 0);
                if (error != "") {
                    cout << error;
                    continue;
                }
                if(match[1]!="D") {
                    cout << "Error: Element " << name << " not found in library\n";
                    continue;
                }if(model=="D" || model=="Z")
                    cout<<controller.addDiode(node1, node2, name, model, circuit);
            } else if (regex_match(input, match, add_ground)) {
                if(match[1]!="GND"){
                    cout<<"Error: Element "<<match[1]<<" not found in library\n";
                    continue;
                }
                cout << controller.addGround(match[2], circuit);
            } else if (regex_match(input, match, delete_ground)) {
                if(match[1]!="GND"){
                    cout<<"Error: Element "<<match[1]<<" not found in library\n";
                    continue;
                }
                cout << controller.deleteGround(match[2], circuit);
            } else if (regex_match(input, match, show_nodes)) {
                controller.showNodes(circuit);
            } else if (regex_match(input, match, list)) {
                controller.list(circuit);
            } else if (regex_match(input, match, list_element)) {
                if(match[1]!="R" && match[1]!="C" && match[1]!="L" && match[1]!="D" && match[1]!="Z") {
                    cout << "Error: Element " << match[1] << " not found in library\n";
                    continue;
                }controller.listElement(match[1], circuit);
            } else if (regex_match(input, match, rename_node)) {
                string old_name = match[1];
                string new_name = match[2];
                string err;
                circuit->renameNode(old_name, new_name, err);
                cout << err;
            } else if (regex_match(input, match, show_details)) {
                controller.showCircuitDetails(circuit);
            } else if (regex_match(input, match, exit)) {
                cout << "Bye Bye!\n";
                return;
            } else
                cout << "Syntax error\n";
        }
    }
};

int main() {
    View view;
    view.run();

    return 0;
}

