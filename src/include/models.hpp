#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.17.0

#include <stan/model/model_header.hpp>

namespace model_idem_namespace {

using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> vector_d;
typedef Eigen::Matrix<double,1,Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> matrix_d;

static int current_statement_begin__;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_idem");
    reader.add_event(102, 102, "end", "model_idem");
    return reader;
}

template <typename T0__, typename T1__, typename T2__>
typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
klpdf(const T0__& err,
          const Eigen::Matrix<T1__, Eigen::Dynamic,1>& res,
          const T2__& h,
          const int& n, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type fun_scalar_t__;
    typedef fun_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        fun_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 6;
        fun_scalar_t__ tmp;
        (void) tmp;  // dummy to suppress unused var warning

        stan::math::initialize(tmp, std::numeric_limits<double>::quiet_NaN());
        stan::math::fill(tmp,DUMMY_VAR__);
        current_statement_begin__ = 7;
        fun_scalar_t__ rst;
        (void) rst;  // dummy to suppress unused var warning

        stan::math::initialize(rst, std::numeric_limits<double>::quiet_NaN());
        stan::math::fill(rst,DUMMY_VAR__);


        current_statement_begin__ = 9;
        stan::math::assign(rst, 0);
        current_statement_begin__ = 10;
        for (int i = 1; i <= n; ++i) {

            current_statement_begin__ = 11;
            stan::math::assign(tmp, ((get_base1(res,i,"res",1) - err) / h));
            current_statement_begin__ = 12;
            stan::math::assign(tmp, pow((tmp / 2),2));
            current_statement_begin__ = 13;
            stan::math::assign(rst, (rst + exp(-(tmp))));
        }
        current_statement_begin__ = 16;
        stan::math::assign(rst, log(((rst / n) / h)));
        current_statement_begin__ = 19;
        return stan::math::promote_scalar<fun_return_scalar_t__>(rst);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct klpdf_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
    operator()(const T0__& err,
          const Eigen::Matrix<T1__, Eigen::Dynamic,1>& res,
          const T2__& h,
          const int& n, std::ostream* pstream__) const {
        return klpdf(err, res, h, n, pstream__);
    }
};

template <bool propto, typename T0__, typename T1__, typename T2__, typename T4__, typename T5__, typename T10__, typename T11__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T4__, typename boost::math::tools::promote_args<T5__, T10__, T11__>::type>::type
cond_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& ymis,
              const Eigen::Matrix<T1__, Eigen::Dynamic,1>& yobs,
              const Eigen::Matrix<T2__, Eigen::Dynamic,Eigen::Dynamic>& coef,
              const int& ny,
              const std::vector<T4__>& mu,
              const std::vector<T5__>& sigma,
              const std::vector<int>& imis,
              const std::vector<int>& inx,
              const int& assumenormal,
              const int& nres,
              const Eigen::Matrix<T10__, Eigen::Dynamic,Eigen::Dynamic>& residual,
              const std::vector<T11__>& h, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T4__, typename boost::math::tools::promote_args<T5__, T10__, T11__>::type>::type fun_scalar_t__;
    typedef fun_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        fun_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 25;
        fun_scalar_t__ rst;
        (void) rst;  // dummy to suppress unused var warning

        stan::math::initialize(rst, std::numeric_limits<double>::quiet_NaN());
        stan::math::fill(rst,DUMMY_VAR__);
        current_statement_begin__ = 26;
        fun_scalar_t__ cmu;
        (void) cmu;  // dummy to suppress unused var warning

        stan::math::initialize(cmu, std::numeric_limits<double>::quiet_NaN());
        stan::math::fill(cmu,DUMMY_VAR__);
        current_statement_begin__ = 27;
        fun_scalar_t__ csigma;
        (void) csigma;  // dummy to suppress unused var warning

        stan::math::initialize(csigma, std::numeric_limits<double>::quiet_NaN());
        stan::math::fill(csigma,DUMMY_VAR__);
        current_statement_begin__ = 28;
        validate_non_negative_index("ally", "ny", ny);
        Eigen::Matrix<fun_scalar_t__,Eigen::Dynamic,1>  ally(static_cast<Eigen::VectorXd::Index>(ny));
        (void) ally;  // dummy to suppress unused var warning

        stan::math::initialize(ally, std::numeric_limits<double>::quiet_NaN());
        stan::math::fill(ally,DUMMY_VAR__);


        current_statement_begin__ = 31;
        for (int i = 1; i <= ny; ++i) {

            current_statement_begin__ = 32;
            if (as_bool(logical_eq(1,get_base1(imis,i,"imis",1)))) {

                current_statement_begin__ = 34;
                stan::math::assign(get_base1_lhs(ally,i,"ally",1), get_base1(ymis,get_base1(inx,i,"inx",1),"ymis",1));
            } else {

                current_statement_begin__ = 37;
                stan::math::assign(get_base1_lhs(ally,i,"ally",1), get_base1(yobs,get_base1(inx,i,"inx",1),"yobs",1));
            }
        }
        current_statement_begin__ = 41;
        stan::math::assign(rst, 0);
        current_statement_begin__ = 42;
        for (int i = 1; i <= ny; ++i) {

            current_statement_begin__ = 43;
            stan::math::assign(csigma, get_base1(sigma,i,"sigma",1));
            current_statement_begin__ = 44;
            stan::math::assign(cmu, get_base1(mu,i,"mu",1));
            current_statement_begin__ = 46;
            if (as_bool(logical_lt(1,i))) {

                current_statement_begin__ = 47;
                stan::math::assign(cmu, (cmu + (get_base1(coef,i,3,"coef",1) * get_base1(ally,(i - 1),"ally",1))));
            }
            current_statement_begin__ = 51;
            if (as_bool(logical_eq(1,assumenormal))) {

                current_statement_begin__ = 52;
                stan::math::assign(rst, (rst + normal_log(get_base1(ally,i,"ally",1),cmu,csigma)));
            } else {

                current_statement_begin__ = 54;
                stan::math::assign(rst, (rst + klpdf((get_base1(ally,i,"ally",1) - cmu),stan::model::rvalue(residual, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(i), stan::model::nil_index_list())), "residual"),get_base1(h,i,"h",1),nres, pstream__)));
            }
        }
        current_statement_begin__ = 57;
        return stan::math::promote_scalar<fun_return_scalar_t__>(rst);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
template <typename T0__, typename T1__, typename T2__, typename T4__, typename T5__, typename T10__, typename T11__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T4__, typename boost::math::tools::promote_args<T5__, T10__, T11__>::type>::type
cond_lpdf(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& ymis,
              const Eigen::Matrix<T1__, Eigen::Dynamic,1>& yobs,
              const Eigen::Matrix<T2__, Eigen::Dynamic,Eigen::Dynamic>& coef,
              const int& ny,
              const std::vector<T4__>& mu,
              const std::vector<T5__>& sigma,
              const std::vector<int>& imis,
              const std::vector<int>& inx,
              const int& assumenormal,
              const int& nres,
              const Eigen::Matrix<T10__, Eigen::Dynamic,Eigen::Dynamic>& residual,
              const std::vector<T11__>& h, std::ostream* pstream__) {
    return cond_lpdf<false>(ymis,yobs,coef,ny,mu,sigma,imis,inx,assumenormal,nres,residual,h, pstream__);
}


struct cond_lpdf_functor__ {
    template <bool propto, typename T0__, typename T1__, typename T2__, typename T4__, typename T5__, typename T10__, typename T11__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T4__, typename boost::math::tools::promote_args<T5__, T10__, T11__>::type>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& ymis,
              const Eigen::Matrix<T1__, Eigen::Dynamic,1>& yobs,
              const Eigen::Matrix<T2__, Eigen::Dynamic,Eigen::Dynamic>& coef,
              const int& ny,
              const std::vector<T4__>& mu,
              const std::vector<T5__>& sigma,
              const std::vector<int>& imis,
              const std::vector<int>& inx,
              const int& assumenormal,
              const int& nres,
              const Eigen::Matrix<T10__, Eigen::Dynamic,Eigen::Dynamic>& residual,
              const std::vector<T11__>& h, std::ostream* pstream__) const {
        return cond_lpdf(ymis, yobs, coef, ny, mu, sigma, imis, inx, assumenormal, nres, residual, h, pstream__);
    }
};

class model_idem : public prob_grad {
private:
    int NY;
    int NOBS;
    vector_d YOBS;
    int NX;
    vector_d X;
    matrix_d COEF;
    vector<int> IMIS;
    vector<int> INX;
    int ASSUMENORMAL;
    int NRES;
    matrix_d RESIDUAL;
    vector<double> H;
    int NMIS;
    vector<double> MU;
    vector<double> SIGMA;
public:
    model_idem(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }

    model_idem(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }

    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning

        current_statement_begin__ = -1;

        static const char* function__ = "model_idem_namespace::model_idem";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        // initialize member variables
        try {
            current_statement_begin__ = 62;
            context__.validate_dims("data initialization", "NY", "int", context__.to_vec());
            NY = int(0);
            vals_i__ = context__.vals_i("NY");
            pos__ = 0;
            NY = vals_i__[pos__++];
            current_statement_begin__ = 63;
            context__.validate_dims("data initialization", "NOBS", "int", context__.to_vec());
            NOBS = int(0);
            vals_i__ = context__.vals_i("NOBS");
            pos__ = 0;
            NOBS = vals_i__[pos__++];
            current_statement_begin__ = 64;
            validate_non_negative_index("YOBS", "(NOBS + 1)", (NOBS + 1));
            context__.validate_dims("data initialization", "YOBS", "vector_d", context__.to_vec((NOBS + 1)));
            validate_non_negative_index("YOBS", "(NOBS + 1)", (NOBS + 1));
            YOBS = vector_d(static_cast<Eigen::VectorXd::Index>((NOBS + 1)));
            vals_r__ = context__.vals_r("YOBS");
            pos__ = 0;
            size_t YOBS_i_vec_lim__ = (NOBS + 1);
            for (size_t i_vec__ = 0; i_vec__ < YOBS_i_vec_lim__; ++i_vec__) {
                YOBS[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 65;
            context__.validate_dims("data initialization", "NX", "int", context__.to_vec());
            NX = int(0);
            vals_i__ = context__.vals_i("NX");
            pos__ = 0;
            NX = vals_i__[pos__++];
            current_statement_begin__ = 66;
            validate_non_negative_index("X", "NX", NX);
            context__.validate_dims("data initialization", "X", "vector_d", context__.to_vec(NX));
            validate_non_negative_index("X", "NX", NX);
            X = vector_d(static_cast<Eigen::VectorXd::Index>(NX));
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_i_vec_lim__ = NX;
            for (size_t i_vec__ = 0; i_vec__ < X_i_vec_lim__; ++i_vec__) {
                X[i_vec__] = vals_r__[pos__++];
            }
            current_statement_begin__ = 67;
            validate_non_negative_index("COEF", "NY", NY);
            validate_non_negative_index("COEF", "(NX + 3)", (NX + 3));
            context__.validate_dims("data initialization", "COEF", "matrix_d", context__.to_vec(NY,(NX + 3)));
            validate_non_negative_index("COEF", "NY", NY);
            validate_non_negative_index("COEF", "(NX + 3)", (NX + 3));
            COEF = matrix_d(static_cast<Eigen::VectorXd::Index>(NY),static_cast<Eigen::VectorXd::Index>((NX + 3)));
            vals_r__ = context__.vals_r("COEF");
            pos__ = 0;
            size_t COEF_m_mat_lim__ = NY;
            size_t COEF_n_mat_lim__ = (NX + 3);
            for (size_t n_mat__ = 0; n_mat__ < COEF_n_mat_lim__; ++n_mat__) {
                for (size_t m_mat__ = 0; m_mat__ < COEF_m_mat_lim__; ++m_mat__) {
                    COEF(m_mat__,n_mat__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 69;
            validate_non_negative_index("IMIS", "NY", NY);
            context__.validate_dims("data initialization", "IMIS", "int", context__.to_vec(NY));
            validate_non_negative_index("IMIS", "NY", NY);
            IMIS = std::vector<int>(NY,int(0));
            vals_i__ = context__.vals_i("IMIS");
            pos__ = 0;
            size_t IMIS_limit_0__ = NY;
            for (size_t i_0__ = 0; i_0__ < IMIS_limit_0__; ++i_0__) {
                IMIS[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 70;
            validate_non_negative_index("INX", "NY", NY);
            context__.validate_dims("data initialization", "INX", "int", context__.to_vec(NY));
            validate_non_negative_index("INX", "NY", NY);
            INX = std::vector<int>(NY,int(0));
            vals_i__ = context__.vals_i("INX");
            pos__ = 0;
            size_t INX_limit_0__ = NY;
            for (size_t i_0__ = 0; i_0__ < INX_limit_0__; ++i_0__) {
                INX[i_0__] = vals_i__[pos__++];
            }
            current_statement_begin__ = 73;
            context__.validate_dims("data initialization", "ASSUMENORMAL", "int", context__.to_vec());
            ASSUMENORMAL = int(0);
            vals_i__ = context__.vals_i("ASSUMENORMAL");
            pos__ = 0;
            ASSUMENORMAL = vals_i__[pos__++];
            current_statement_begin__ = 74;
            context__.validate_dims("data initialization", "NRES", "int", context__.to_vec());
            NRES = int(0);
            vals_i__ = context__.vals_i("NRES");
            pos__ = 0;
            NRES = vals_i__[pos__++];
            current_statement_begin__ = 75;
            validate_non_negative_index("RESIDUAL", "NRES", NRES);
            validate_non_negative_index("RESIDUAL", "NY", NY);
            context__.validate_dims("data initialization", "RESIDUAL", "matrix_d", context__.to_vec(NRES,NY));
            validate_non_negative_index("RESIDUAL", "NRES", NRES);
            validate_non_negative_index("RESIDUAL", "NY", NY);
            RESIDUAL = matrix_d(static_cast<Eigen::VectorXd::Index>(NRES),static_cast<Eigen::VectorXd::Index>(NY));
            vals_r__ = context__.vals_r("RESIDUAL");
            pos__ = 0;
            size_t RESIDUAL_m_mat_lim__ = NRES;
            size_t RESIDUAL_n_mat_lim__ = NY;
            for (size_t n_mat__ = 0; n_mat__ < RESIDUAL_n_mat_lim__; ++n_mat__) {
                for (size_t m_mat__ = 0; m_mat__ < RESIDUAL_m_mat_lim__; ++m_mat__) {
                    RESIDUAL(m_mat__,n_mat__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 76;
            validate_non_negative_index("H", "NY", NY);
            context__.validate_dims("data initialization", "H", "double", context__.to_vec(NY));
            validate_non_negative_index("H", "NY", NY);
            H = std::vector<double>(NY,double(0));
            vals_r__ = context__.vals_r("H");
            pos__ = 0;
            size_t H_limit_0__ = NY;
            for (size_t i_0__ = 0; i_0__ < H_limit_0__; ++i_0__) {
                H[i_0__] = vals_r__[pos__++];
            }

            // validate, data variables
            current_statement_begin__ = 62;
            check_greater_or_equal(function__,"NY",NY,1);
            current_statement_begin__ = 63;
            check_greater_or_equal(function__,"NOBS",NOBS,0);
            current_statement_begin__ = 64;
            current_statement_begin__ = 65;
            check_greater_or_equal(function__,"NX",NX,1);
            current_statement_begin__ = 66;
            current_statement_begin__ = 67;
            current_statement_begin__ = 69;
            for (int k0__ = 0; k0__ < NY; ++k0__) {
                check_greater_or_equal(function__,"IMIS[k0__]",IMIS[k0__],0);
                check_less_or_equal(function__,"IMIS[k0__]",IMIS[k0__],1);
            }
            current_statement_begin__ = 70;
            for (int k0__ = 0; k0__ < NY; ++k0__) {
                check_greater_or_equal(function__,"INX[k0__]",INX[k0__],1);
            }
            current_statement_begin__ = 73;
            check_greater_or_equal(function__,"ASSUMENORMAL",ASSUMENORMAL,0);
            check_less_or_equal(function__,"ASSUMENORMAL",ASSUMENORMAL,1);
            current_statement_begin__ = 74;
            check_greater_or_equal(function__,"NRES",NRES,1);
            current_statement_begin__ = 75;
            current_statement_begin__ = 76;
            // initialize data variables
            current_statement_begin__ = 81;
            NMIS = int(0);
            stan::math::fill(NMIS, std::numeric_limits<int>::min());
            current_statement_begin__ = 82;
            validate_non_negative_index("MU", "NY", NY);
            MU = std::vector<double>(NY,double(0));
            stan::math::fill(MU,DUMMY_VAR__);
            current_statement_begin__ = 83;
            validate_non_negative_index("SIGMA", "NY", NY);
            SIGMA = std::vector<double>(NY,double(0));
            stan::math::fill(SIGMA,DUMMY_VAR__);

            current_statement_begin__ = 85;
            stan::math::assign(NMIS, (NY - NOBS));
            current_statement_begin__ = 86;
            for (int i = 1; i <= NY; ++i) {

                current_statement_begin__ = 87;
                stan::math::assign(get_base1_lhs(SIGMA,i,"SIGMA",1), get_base1(COEF,i,1,"COEF",1));
                current_statement_begin__ = 88;
                stan::math::assign(get_base1_lhs(MU,i,"MU",1), get_base1(COEF,i,2,"COEF",1));
                current_statement_begin__ = 90;
                for (int k = 1; k <= NX; ++k) {

                    current_statement_begin__ = 91;
                    stan::math::assign(get_base1_lhs(MU,i,"MU",1), (get_base1(MU,i,"MU",1) + (get_base1(X,k,"X",1) * get_base1(COEF,i,(3 + k),"COEF",1))));
                }
            }

            // validate transformed data
            current_statement_begin__ = 81;
            check_greater_or_equal(function__,"NMIS",NMIS,0);
            current_statement_begin__ = 82;
            current_statement_begin__ = 83;

            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 97;
            validate_non_negative_index("YMIS", "NMIS", NMIS);
            num_params_r__ += NMIS;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    ~model_idem() { }


    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        stan::io::writer<double> writer__(params_r__,params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;

        if (!(context__.contains_r("YMIS")))
            throw std::runtime_error("variable YMIS missing");
        vals_r__ = context__.vals_r("YMIS");
        pos__ = 0U;
        validate_non_negative_index("YMIS", "NMIS", NMIS);
        context__.validate_dims("initialization", "YMIS", "vector_d", context__.to_vec(NMIS));
        vector_d YMIS(static_cast<Eigen::VectorXd::Index>(NMIS));
        for (int j1__ = 0U; j1__ < NMIS; ++j1__)
            YMIS(j1__) = vals_r__[pos__++];
        try {
            writer__.vector_unconstrain(YMIS);
        } catch (const std::exception& e) { 
            throw std::runtime_error(std::string("Error transforming variable YMIS: ") + e.what());
        }

        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }


    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(vector<T__>& params_r__,
                 vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {

        T__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;

        try {
            // model parameters
            stan::io::reader<T__> in__(params_r__,params_i__);

            Eigen::Matrix<T__,Eigen::Dynamic,1>  YMIS;
            (void) YMIS;  // dummy to suppress unused var warning
            if (jacobian__)
                YMIS = in__.vector_constrain(NMIS,lp__);
            else
                YMIS = in__.vector_constrain(NMIS);


            // transformed parameters



            // validate transformed parameters

            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning

            // model body

            current_statement_begin__ = 101;
            lp_accum__.add(cond_lpdf<propto__>(YMIS, YOBS, COEF, NY, MU, SIGMA, IMIS, INX, ASSUMENORMAL, NRES, RESIDUAL, H, pstream__));

        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }

        lp_accum__.add(lp__);
        return lp_accum__.sum();

    } // log_prob()

    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }


    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("YMIS");
    }


    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(NMIS);
        dimss__.push_back(dims__);
    }

    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        vars__.resize(0);
        stan::io::reader<double> in__(params_r__,params_i__);
        static const char* function__ = "model_idem_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        vector_d YMIS = in__.vector_constrain(NMIS);
            for (int k_0__ = 0; k_0__ < NMIS; ++k_0__) {
            vars__.push_back(YMIS[k_0__]);
            }

        if (!include_tparams__) return;
        // declare and define transformed parameters
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;

        double DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

        try {



            // validate transformed parameters

            // write transformed parameters

            if (!include_gqs__) return;
            // declare and define generated quantities



            // validate generated quantities

            // write generated quantities
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }

    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng,params_r_vec,params_i_vec,vars_vec,include_tparams,include_gqs,pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    static std::string model_name() {
        return "model_idem";
    }


    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= NMIS; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "YMIS" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }


    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        for (int k_0__ = 1; k_0__ <= NMIS; ++k_0__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "YMIS" << '.' << k_0__;
            param_names__.push_back(param_name_stream__.str());
        }

        if (!include_gqs__ && !include_tparams__) return;

        if (!include_gqs__) return;
    }

}; // model

}




#endif
