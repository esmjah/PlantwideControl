within networkModels;
package networkComponents 
  replaceable function maxFunc
    input Real x1;
    input Real x2;
    output Real y;
  algorithm
    y := max(x1, x2);
  end maxFunc;


  replaceable function minFunc
    input Real x1;
    input Real x2;
    output Real y;
  algorithm
    y := min(x1, x2);
  end minFunc;
end networkComponents;
