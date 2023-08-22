(* Content-type: application/vnd.wolfram.mathematica *)

(*** Wolfram Notebook File ***)
(* http://www.wolfram.com/nb *)

(* CreatedBy='Mathematica 13.2' *)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[       158,          7]
NotebookDataLength[     11088,        307]
NotebookOptionsPosition[      7415,        230]
NotebookOutlinePosition[      9068,        268]
CellTagsIndexPosition[      8987,        263]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell[BoxData[
 RowBox[{"Integrate", "[", 
  RowBox[{
   RowBox[{"1", "+", 
    RowBox[{"x", "^", "2"}], " ", "+", " ", 
    RowBox[{"y", "^", "2"}], " ", "-", 
    RowBox[{"2", "*", 
     RowBox[{"y", "^", "3"}]}]}], ",", 
   RowBox[{"{", 
    RowBox[{"x", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"y", ",", "0", ",", "1"}], "}"}]}], "]"}]], "Input",
 CellChangeTimes->{{3.887622708837659*^9, 3.88762282725694*^9}},
 CellTags->"Integrate_templates",
 CellLabel->"In[2]:=",ExpressionUUID->"deb446b3-915b-4e83-a6b5-b6d6c5ac36f9"],

Cell[BoxData[
 FractionBox["7", "6"]], "Output",
 CellChangeTimes->{{3.8876228113235292`*^9, 3.887622828931079*^9}},
 CellTags->"Integrate_templates",
 CellLabel->"Out[2]=",ExpressionUUID->"d24543a3-bbbc-42e4-88b6-7bcea29e2c85"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Integrate", "[", 
  RowBox[{
   RowBox[{"1", "/", "x"}], ",", 
   RowBox[{"{", 
    RowBox[{"x", ",", "0", ",", "1"}], "}"}]}], "]"}]], "Input",
 CellChangeTimes->{{3.8876230617831755`*^9, 3.8876230881464753`*^9}},
 CellTags->"Integrate_templates",
 CellLabel->"In[3]:=",ExpressionUUID->"398d4a3c-677e-4cbf-a50f-d9c523050fd7"],

Cell[BoxData[
 TemplateBox[{
  "Integrate", "idiv", 
   "\"Integral of \\!\\(\\*FractionBox[\\\"1\\\", \\\"x\\\"]\\) does not \
converge on \\!\\(\\*RowBox[{\\\"{\\\", RowBox[{\\\"0\\\", \\\",\\\", \\\"1\\\
\"}], \\\"}\\\"}]\\).\"", 2, 3, 1, 27705338991265594300, "Local"},
  "MessageTemplate"]], "Message", "MSG",
 CellChangeTimes->{3.8876230906597548`*^9},
 CellTags->"Integrate_templates",
 CellLabel->
  "During evaluation of \
In[3]:=",ExpressionUUID->"10f2987a-8314-4719-bb82-e70baeeb8638"],

Cell[BoxData[
 RowBox[{
  SubsuperscriptBox["\[Integral]", "0", "1"], 
  RowBox[{
   FractionBox["1", "x"], 
   RowBox[{"\[DifferentialD]", "x"}]}]}]], "Output",
 CellChangeTimes->{3.88762309068479*^9},
 CellTags->"Integrate_templates",
 CellLabel->"Out[3]=",ExpressionUUID->"4c25ba4c-c8ff-4466-8276-f1633088f04f"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"N", "[", 
  RowBox[{"Integrate", "[", 
   RowBox[{
    RowBox[{"1", "/", "x"}], ",", 
    RowBox[{"{", 
     RowBox[{"x", ",", 
      SuperscriptBox["10", 
       RowBox[{"-", "300"}]], ",", "1"}], "}"}]}], "]"}], "]"}]], "Input",
 CellChangeTimes->{{3.8876230617831755`*^9, 3.8876230881464753`*^9}, {
  3.887623407513091*^9, 3.88762351476097*^9}, {3.887623598025052*^9, 
  3.8876238215608864`*^9}, {3.8876346654527607`*^9, 3.887634713886326*^9}},
 CellTags->"Integrate_templates",
 CellLabel->"In[46]:=",ExpressionUUID->"fbbbc51b-b243-4552-86fa-84633eddb0bd"],

Cell[BoxData["690.7755278982137`"], "Output",
 CellChangeTimes->{
  3.88762309068479*^9, {3.8876234112540984`*^9, 3.8876235169726915`*^9}, {
   3.8876236099743185`*^9, 3.8876236166250257`*^9}, {3.887623671295089*^9, 
   3.8876238233718805`*^9}, {3.887634674008326*^9, 3.887634714965717*^9}},
 CellTags->"Integrate_templates",
 CellLabel->"Out[46]=",ExpressionUUID->"a10830a0-1cbd-46cb-86e0-c9e97408a042"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"N", "[", 
  RowBox[{"Integrate", "[", 
   RowBox[{
    RowBox[{"x", "*", "y"}], ",", 
    RowBox[{"{", 
     RowBox[{"y", ",", 
      RowBox[{"-", "3"}], ",", "0"}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"x", ",", "2", ",", 
      RowBox[{"2", "-", "y"}]}], "}"}]}], "]"}], "]"}]], "Input",
 CellChangeTimes->{{3.88762458502017*^9, 3.8876246421920543`*^9}},
 CellTags->"Integrate_templates",
 CellLabel->"In[35]:=",ExpressionUUID->"009c8805-78b5-4f94-9eb0-684c114f099b"],

Cell[BoxData[
 RowBox[{"-", "28.125`"}]], "Output",
 CellChangeTimes->{{3.887624620122362*^9, 3.8876246434417953`*^9}},
 CellTags->"Integrate_templates",
 CellLabel->"Out[35]=",ExpressionUUID->"7bdd1a6a-720f-4bcb-adff-daaf3ee72074"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Integrate", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1", "-", "x"}], ")"}], "*", 
    RowBox[{"(", 
     RowBox[{"1", "-", "y"}], ")"}], "*", "x", "*", 
    RowBox[{"(", 
     RowBox[{"1", "-", "y"}], ")"}]}], ",", 
   RowBox[{"{", 
    RowBox[{"x", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"y", ",", "0", ",", "1"}], "}"}]}], "]"}]], "Input",
 CellChangeTimes->{{3.887624795203742*^9, 3.887624819032911*^9}, {
  3.8876249688028245`*^9, 3.887624995623569*^9}},
 CellLabel->"In[37]:=",ExpressionUUID->"e3c51503-25bb-41af-aa8b-cb68b1ac4e33"],

Cell[BoxData[
 FractionBox["1", "18"]], "Output",
 CellChangeTimes->{3.8876248214262514`*^9, 3.887624997223088*^9},
 CellLabel->"Out[37]=",ExpressionUUID->"dcabcde3-f4e5-4559-83ee-cf7f2c6c370b"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Integrate", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1", "-", 
      RowBox[{"2", "x"}]}], ")"}], 
    RowBox[{"(", 
     RowBox[{"1", "+", 
      SuperscriptBox["y", "2"], "-", 
      RowBox[{"2", "y"}]}], ")"}]}], ",", 
   RowBox[{"{", 
    RowBox[{"x", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"y", ",", "0", ",", "1"}], "}"}]}], "]"}]], "Input",
 CellChangeTimes->{{3.8876253772430286`*^9, 3.887625434831303*^9}},
 CellLabel->"In[38]:=",ExpressionUUID->"c5e84e59-5c25-4953-8812-881efae94422"],

Cell[BoxData["0"], "Output",
 CellChangeTimes->{3.8876254397834806`*^9},
 CellLabel->"Out[38]=",ExpressionUUID->"8e97f90c-d7c9-45e9-bcad-3dc389ffec08"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Integrate", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"(", 
     RowBox[{"x", "-", 
      SuperscriptBox["x", "2"]}], ")"}], 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"2", "y"}], "-", "2"}], ")"}]}], ",", 
   RowBox[{"{", 
    RowBox[{"x", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"y", ",", "0", ",", "1"}], "}"}]}], "]"}]], "Input",
 CellChangeTimes->{{3.8876254898679476`*^9, 3.8876254991117606`*^9}},
 CellLabel->"In[39]:=",ExpressionUUID->"b57dd3d9-b4f7-4387-8215-2a5915a6cc53"],

Cell[BoxData[
 RowBox[{"-", 
  FractionBox["1", "6"]}]], "Output",
 CellChangeTimes->{3.887625501573184*^9},
 CellLabel->"Out[39]=",ExpressionUUID->"f3b6b4cb-0b5e-453f-8eb9-23f0a88e64e5"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"N", "[", 
  RowBox[{"Integrate", "[", 
   RowBox[{
    RowBox[{"x", "*", "y"}], ",", 
    RowBox[{"{", 
     RowBox[{"x", ",", "2", ",", "5"}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"y", ",", 
      RowBox[{"-", "3"}], ",", 
      RowBox[{"2", "-", "x"}]}], "}"}]}], "]"}], "]"}]], "Input",
 CellChangeTimes->{{3.8876258676100335`*^9, 3.8876258996473227`*^9}},
 CellLabel->"In[40]:=",ExpressionUUID->"a91c1219-515e-4873-864d-e7283a5e0f07"],

Cell[BoxData[
 RowBox[{"-", "28.125`"}]], "Output",
 CellChangeTimes->{3.887625903314657*^9},
 CellLabel->"Out[40]=",ExpressionUUID->"f83ccefc-74dc-4f3e-bb08-dc4ce9fe3673"]
}, Open  ]]
},
WindowSize->{1440., 837.75},
WindowMargins->{{-6, Automatic}, {Automatic, -6}},
FrontEndVersion->"13.2 for Microsoft Windows (64-bit) (January 30, 2023)",
StyleDefinitions->"Default.nb",
ExpressionUUID->"2ac4ebdc-172f-4068-a066-e4ed696139c7"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{
 "Integrate_templates"->{
  Cell[580, 22, 550, 14, 28, "Input",ExpressionUUID->"deb446b3-915b-4e83-a6b5-b6d6c5ac36f9",
   CellTags->"Integrate_templates"],
  Cell[1133, 38, 228, 4, 48, "Output",ExpressionUUID->"d24543a3-bbbc-42e4-88b6-7bcea29e2c85",
   CellTags->"Integrate_templates"],
  Cell[1398, 47, 350, 8, 28, "Input",ExpressionUUID->"398d4a3c-677e-4cbf-a50f-d9c523050fd7",
   CellTags->"Integrate_templates"],
  Cell[1751, 57, 496, 11, 38, "Message",ExpressionUUID->"10f2987a-8314-4719-bb82-e70baeeb8638",
   CellTags->"Integrate_templates"],
  Cell[2250, 70, 314, 8, 49, "Output",ExpressionUUID->"4c25ba4c-c8ff-4466-8276-f1633088f04f",
   CellTags->"Integrate_templates"],
  Cell[2601, 83, 584, 13, 31, "Input",ExpressionUUID->"fbbbc51b-b243-4552-86fa-84633eddb0bd",
   CellTags->"Integrate_templates"],
  Cell[3188, 98, 404, 6, 32, "Output",ExpressionUUID->"a10830a0-1cbd-46cb-86e0-c9e97408a042",
   CellTags->"Integrate_templates"],
  Cell[3629, 109, 500, 13, 28, "Input",ExpressionUUID->"009c8805-78b5-4f94-9eb0-684c114f099b",
   CellTags->"Integrate_templates"],
  Cell[4132, 124, 232, 4, 32, "Output",ExpressionUUID->"7bdd1a6a-720f-4bcb-adff-daaf3ee72074",
   CellTags->"Integrate_templates"]}
 }
*)
(*CellTagsIndex
CellTagsIndex->{
 {"Integrate_templates", 7782, 241}
 }
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[580, 22, 550, 14, 28, "Input",ExpressionUUID->"deb446b3-915b-4e83-a6b5-b6d6c5ac36f9",
 CellTags->"Integrate_templates"],
Cell[1133, 38, 228, 4, 48, "Output",ExpressionUUID->"d24543a3-bbbc-42e4-88b6-7bcea29e2c85",
 CellTags->"Integrate_templates"]
}, Open  ]],
Cell[CellGroupData[{
Cell[1398, 47, 350, 8, 28, "Input",ExpressionUUID->"398d4a3c-677e-4cbf-a50f-d9c523050fd7",
 CellTags->"Integrate_templates"],
Cell[1751, 57, 496, 11, 38, "Message",ExpressionUUID->"10f2987a-8314-4719-bb82-e70baeeb8638",
 CellTags->"Integrate_templates"],
Cell[2250, 70, 314, 8, 49, "Output",ExpressionUUID->"4c25ba4c-c8ff-4466-8276-f1633088f04f",
 CellTags->"Integrate_templates"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2601, 83, 584, 13, 31, "Input",ExpressionUUID->"fbbbc51b-b243-4552-86fa-84633eddb0bd",
 CellTags->"Integrate_templates"],
Cell[3188, 98, 404, 6, 32, "Output",ExpressionUUID->"a10830a0-1cbd-46cb-86e0-c9e97408a042",
 CellTags->"Integrate_templates"]
}, Open  ]],
Cell[CellGroupData[{
Cell[3629, 109, 500, 13, 28, "Input",ExpressionUUID->"009c8805-78b5-4f94-9eb0-684c114f099b",
 CellTags->"Integrate_templates"],
Cell[4132, 124, 232, 4, 32, "Output",ExpressionUUID->"7bdd1a6a-720f-4bcb-adff-daaf3ee72074",
 CellTags->"Integrate_templates"]
}, Open  ]],
Cell[CellGroupData[{
Cell[4401, 133, 606, 16, 28, "Input",ExpressionUUID->"e3c51503-25bb-41af-aa8b-cb68b1ac4e33"],
Cell[5010, 151, 194, 3, 48, "Output",ExpressionUUID->"dcabcde3-f4e5-4559-83ee-cf7f2c6c370b"]
}, Open  ]],
Cell[CellGroupData[{
Cell[5241, 159, 563, 16, 31, "Input",ExpressionUUID->"c5e84e59-5c25-4953-8812-881efae94422"],
Cell[5807, 177, 151, 2, 32, "Output",ExpressionUUID->"8e97f90c-d7c9-45e9-bcad-3dc389ffec08"]
}, Open  ]],
Cell[CellGroupData[{
Cell[5995, 184, 533, 15, 31, "Input",ExpressionUUID->"b57dd3d9-b4f7-4387-8215-2a5915a6cc53"],
Cell[6531, 201, 187, 4, 48, "Output",ExpressionUUID->"f3b6b4cb-0b5e-453f-8eb9-23f0a88e64e5"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6755, 210, 469, 12, 28, "Input",ExpressionUUID->"a91c1219-515e-4873-864d-e7283a5e0f07"],
Cell[7227, 224, 172, 3, 32, "Output",ExpressionUUID->"f83ccefc-74dc-4f3e-bb08-dc4ce9fe3673"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

