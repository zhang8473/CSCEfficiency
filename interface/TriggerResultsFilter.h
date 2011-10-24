#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"

// forward declaration
namespace triggerExpression {
  class Evaluator;
}

//
// class declaration
//

class TriggerResultsFilter {
public:
  explicit TriggerResultsFilter(const edm::ParameterSet &);
  ~TriggerResultsFilter();
  virtual bool TriggerFilter(const edm::Event &, const edm::EventSetup &);

private:
  /// parse the logical expression into functionals
  void parse(const std::string & expression);
  void parse(const std::vector<std::string> & expressions);

  /// evaluator for the trigger condition
  triggerExpression::Evaluator * m_expression;

  /// cache some data from the Event for faster access by the m_expression
  triggerExpression::Data m_eventCache;
};

//
// constructors and destructor
//
TriggerResultsFilter::TriggerResultsFilter(const edm::ParameterSet & config) :
  m_expression(0),
  m_eventCache(config)
{
  const std::vector<std::string> & expressions = config.getParameter<std::vector<std::string> >("triggerConditions");
  parse( expressions );
}

TriggerResultsFilter::~TriggerResultsFilter()
{
  delete m_expression;
}

void TriggerResultsFilter::parse(const std::vector<std::string> & expressions) {
  // parse the logical expressions into functionals
  if (expressions.size() == 0) {
    edm::LogWarning("Configuration") << "Empty trigger results expression";
  } else if (expressions.size() == 1) {
    parse( expressions[0] );
  } else {
    std::stringstream expression;
    expression << "(" << expressions[0] << ")";
    for (unsigned int i = 1; i < expressions.size(); ++i)
      expression << " OR (" << expressions[i] << ")";
    parse( expression.str() );
  }
}

void TriggerResultsFilter::parse(const std::string & expression) {
  // parse the logical expressions into functionals
  m_expression = triggerExpression::parse( expression );

  // check if the expressions were parsed correctly
  if (not m_expression)
    edm::LogWarning("Configuration") << "Couldn't parse trigger results expression \"" << expression << "\"";
}

bool TriggerResultsFilter::TriggerFilter(const edm::Event & event, const edm::EventSetup & setup)
{
  if (not m_expression)
    // no valid expression has been parsed
    return false;

  if (not m_eventCache.setEvent(event, setup))
    // couldn't properly access all information from the Event
    return false;

  // if the L1 or HLT configurations have changed, (re)initialize the filters (including during the first event)
  if (m_eventCache.configurationUpdated()) {
    m_expression->init(m_eventCache);

    // log the expanded configuration
    edm::LogInfo("Configuration") << "TriggerResultsFilter configuration updated: " << *m_expression;
  }

  // run the trigger results filter
  return (*m_expression)(m_eventCache);
}
