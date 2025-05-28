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
    vector <Element*> getConnectedElements() {return connected_elements;}
    string getName(){return name;}

    void addConnectedElement(Element *e){
        connected_elements.push_back(e);
    }

    void removeConnectedElement(Element* e){
        auto x=find(connected_elements.begin(), connected_elements.end(), e);
        if(x!=connected_elements.end())
            connected_elements.erase(x);
    }
};

class Circuit {
protected:
    vector<Element *> elements;
    vector<Node *> nodes;
    map<string, Node *> node_access;
public:
    vector<Element *> &getElements() {return elements;}
    map<string, Node*> getNodeAccess() {return node_access;}

    Node* getCreateNode(string name){
        if(!node_access.count(name)){
            Node* newnode= new Node(name);
            nodes.push_back(newnode);
            node_access[name]=newnode;
        }
        return node_access[name];
    }
    void addElement(Element* element){
        elements.push_back(element);
    }
    void removeElement(Element* element){
        auto x= find(elements.begin(), elements.end(), element);
        if(x != elements.end())
            elements.erase(x);
    }
};

class Element {
protected:
    string name;
    double value;
    Node *node1;
    Node *node2;
public:
    Element(string name_, double value_, Node* n1, Node* n2):
            name(name_), value(value_), node1(n1), node2(n2) {}

    string getName() {
        return name;
    }
    virtual string getType() = 0;
    double getValue()  {return value;}
    Node* getFirstNode() {return node1;}
    Node* getSecondNode() {return node2;}
    virtual ~Element() = default;
};

class Resistor : public Element{
public:
    Resistor(string name_, double value_, Node* n1, Node* n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override {return "Resistor";}
};

class Capacitor : public Element{
public:
    Capacitor(string name_, double value_, Node* n1, Node* n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override {return "Capacitor";}
};

class Inductor : public Element{
public:
    Inductor(string name_, double value_, Node* n1, Node* n2)
            : Element(name_, value_, n1, n2) {}

    string getType() override {return "Inductor";}
};

class Controller {
public:
    string handleError(string name, string number, Circuit *circuit, double value) {
        if (value <= 0)
            return "Error: Resistance cannot be zero or negative\n";
        for (auto element: circuit->getElements()) {
            if(!element)
                continue;
            if (element->getName() == name) {
                switch (name[0]) {
                    case 'R':
                        return "Resistor " + name + " already exists in the circuit\n";
                    case 'C':
                        return "Capacitor " + name + " already exists in the circuit\n";
                    case 'L':
                        return "Inductor " + name + " already exists in the circuit\n";
                }
            }
        }
        return "";
    }
    string addNewElement(string node1, string node2, string name, double value, Circuit* circuit){
        Node* n1=circuit->getCreateNode(node1);
        Node* n2=circuit->getCreateNode(node2);
        Element* element= nullptr;
        switch(name[0]){
            case 'R':
                element=new Resistor(name, value, n1, n2);
                break;
            case 'L':
                element=new Inductor(name, value, n1, n2);
                break;
            case 'C':
                element=new Capacitor(name, value, n1, n2);
                break;
        }
        circuit->addElement(element);
        element->getFirstNode()->addConnectedElement(element);
        element->getSecondNode()->addConnectedElement(element);
        return element->getType() + " added successfully!\n";
    }
    string removeElement(string name, Circuit* circuit){
        auto& elements=circuit->getElements();
        for(auto it= elements.begin(); it!=elements.end(); it++){
            if((*it)->getName() == name){
                (*it)->getFirstNode()->removeConnectedElement(*it);
                (*it)->getSecondNode()->removeConnectedElement(*it);
                Element* deleting=*it;
                elements.erase(it);
                delete deleting;
                return name + " removed successfully!\n";
            }
        }
        switch (name[0]) {
            case'R':
                return "Error: Cannot delete resistor; component not found\n";
            case 'L':
                return "Error: Cannot delete inductor; component not found\n";
            case 'C':
                return "Error: Cannot delete capacitor; component not found\n";
        }
    }
    void showCircuitDetails(Circuit* circuit){
        cout<<"Circuit Details:\nElements:\n";
        for(auto element: circuit->getElements()){
            cout<<"type: "<<element->getType()<<"   name: "<<element->getName()<<"  value: "<<element->getValue()
                <<"  node1: "<<element->getFirstNode()->getName()<<"  node2: "<<element->getSecondNode()->getName()<<endl;
        }cout<<"Node Details:\n";
        for(auto node: circuit->getNodeAccess()){
            cout<<node.first<<" : "<<endl;
            if(node.second->getConnectedElements().size()==0)
                cout<<"No element is connected to this node!\n";
            else{
                for(auto element: node.second->getConnectedElements())
                    cout<<element->getName()<<endl;
            }
        }
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
        regex show_details(R"(^\s*show\s+details\s*$)");
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
                    cout<< "Error: Invalid numeric value\n";
                    continue;
                }
                string error=controller.handleError(name, number, circuit, value);
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
                cout<<controller.addNewElement(node1, node2, name, value, circuit);
                cout << name << " " << node1 << " " << node2 << " " << number << " " << unit << " "<<value<< "\n";

            } else if(regex_match(input, match, remove_element)){
                if (match[1] != "R" && match[1] != "L" && match[1] != "C") {
                    cout << "Element " << match[1] << " not found in library\n";
                    continue;
                }
                cout<<controller.removeElement(match[1].str() + match[2].str(), circuit);
            } else if(regex_match(input, match, show_details)){
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