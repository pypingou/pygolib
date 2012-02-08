%{!?python_sitelib: %global python_sitelib %(%{__python} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())")}

Name:           pygolib
Version:        0.1.0
Release:        1%{?dist}
Summary:        A python library to manipulate GO graph

License:        BSD
URL:            https://github.com/PBR/pygolib
Source0:        pygolib-0.1.0.tar.gz

BuildArch:      noarch
BuildRequires:  python-devel

%description
A python library to manipulate the Gene Ontology graph.

This library offer you a way to download, load the GO graph.
Distance between two GO terms can be computed taking into
account both the distance in the tree and the distance to
last common ancester.

%prep
%setup -q


%build
%{__python} setup.py build


%install
rm -rf $RPM_BUILD_ROOT
%{__python} setup.py install -O1 --skip-build --root $RPM_BUILD_ROOT

 
%files
%doc COPYING README
# For noarch packages: sitelib
%{python_sitelib}/*
%{_bindir}/goutil

%changelog
* Wed Feb 08 2012 Pierre-Yves Chibon <pingou@pingoured.fr> - 0.1.0-1
- Initial spec file according to Fedora's guidelines
